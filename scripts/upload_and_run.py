#!/usr/bin/env python3
"""
Upload pipeline_data + scripts to VM and launch the scoring pipeline.

Usage:
    python upload_and_run.py [--vm_host IP] [--dry_run]

Default VM host: 10.110.5.102 (update if IP changed)
"""
import argparse
import io
import os
import sys
import time
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

try:
    import paramiko
    from paramiko import SFTPClient
except ImportError:
    print("ERROR: pip install paramiko")
    sys.exit(1)

LOCAL_BASE  = Path(r"G:\My Drive\PEN-STACK\PAPER_4\pen-assemble")
KEY_FILE    = r"C:\Users\ANEES AHMED\.ssh\google_compute_engine"
USERNAME    = "anees_22phd0670"

REMOTE_DATA    = "/home/anees_22phd0670/pen_pipeline_data"
REMOTE_SCRIPTS = "/home/anees_22phd0670/pen_pipeline_scripts"
REMOTE_RESULTS = "/home/anees_22phd0670/pen_pipeline_results"
REMOTE_LOG     = "/home/anees_22phd0670/pen_pipeline_results/logs/pipeline_run.log"


def try_connect(host: str, timeout: int = 10) -> paramiko.SSHClient:
    """Try RSA key first, then Ed25519."""
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    for key_type in [paramiko.RSAKey, paramiko.Ed25519Key, paramiko.ECDSAKey]:
        try:
            key = key_type.from_private_key_file(KEY_FILE)
            client.connect(host, username=USERNAME, pkey=key, timeout=timeout)
            print(f"  Connected to {host} (key type: {key_type.__name__})")
            return client
        except Exception:
            continue

    # Try without explicit key type
    try:
        client.connect(host, username=USERNAME, key_filename=KEY_FILE, timeout=timeout)
        print(f"  Connected to {host}")
        return client
    except Exception as e:
        raise ConnectionError(f"Cannot connect to {host}: {e}")


def exec_cmd(client: paramiko.SSHClient, cmd: str, print_output: bool = True) -> tuple[int, str]:
    _, stdout, stderr = client.exec_command(cmd)
    out = stdout.read().decode('utf-8', 'replace').strip()
    err = stderr.read().decode('utf-8', 'replace').strip()
    rc  = stdout.channel.recv_exit_status()
    if print_output and out:
        print(out)
    if err:
        print(f"  STDERR: {err[:500]}")
    return rc, out


def sftp_upload_dir(sftp: SFTPClient, local_dir: Path, remote_dir: str):
    """Recursively upload local_dir to remote_dir."""
    try:
        sftp.mkdir(remote_dir)
    except Exception:
        pass
    for item in sorted(local_dir.iterdir()):
        remote_path = f"{remote_dir}/{item.name}"
        if item.is_dir():
            sftp_upload_dir(sftp, item, remote_path)
        else:
            size_kb = item.stat().st_size // 1024
            print(f"  Uploading {item.relative_to(LOCAL_BASE)} ({size_kb} KB)...")
            sftp.put(str(item), remote_path)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--vm_host",  default="10.110.5.102", help="VM IP address")
    p.add_argument("--dry_run",  action="store_true")
    p.add_argument("--skip_upload", action="store_true", help="Skip upload, just run")
    args = p.parse_args()

    print(f"\n{'='*60}")
    print("PEN-ASSEMBLE Pipeline Upload & Run")
    print(f"{'='*60}")
    print(f"VM host:  {args.vm_host}")
    print(f"Dry run:  {args.dry_run}")

    # ── Connect ────────────────────────────────────────────────────────────────
    print(f"\nConnecting to {args.vm_host}...")
    client = try_connect(args.vm_host)
    rc, out = exec_cmd(client, "hostname && uname -r && nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1")
    print(f"  VM: {out}")

    if args.dry_run:
        print("\n[DRY RUN] Would upload and execute. Exiting.")
        client.close()
        return

    sftp = client.open_sftp()

    # ── Upload pipeline_data ───────────────────────────────────────────────────
    if not args.skip_upload:
        print(f"\nUploading pipeline_data to {REMOTE_DATA}...")
        sftp_upload_dir(sftp, LOCAL_BASE / "pipeline_data", REMOTE_DATA)

        print(f"\nUploading pipeline scripts to {REMOTE_SCRIPTS}...")
        try:
            sftp.mkdir(REMOTE_SCRIPTS)
        except Exception:
            pass
        script = LOCAL_BASE / "scripts/run_pipeline_direct.py"
        sftp.put(str(script), f"{REMOTE_SCRIPTS}/run_pipeline_direct.py")
        print(f"  Uploaded run_pipeline_direct.py")

    # ── Verify upload ──────────────────────────────────────────────────────────
    print("\nVerifying uploads...")
    exec_cmd(client, f"ls {REMOTE_DATA}/designs/ && ls {REMOTE_DATA}/structures/")

    # ── Check dependencies ─────────────────────────────────────────────────────
    print("\nChecking VM dependencies...")
    exec_cmd(client, "python3 -c \"import pandas, pyarrow; print('pandas/pyarrow OK')\"")
    exec_cmd(client, "python3 -c \"import mech_class; print('mech_class:', mech_class.__version__)\" 2>&1")
    exec_cmd(client, "python3 -c \"import pen_score; print('pen_score:', pen_score.__version__)\" 2>&1")
    exec_cmd(client, "python3 -c \"import pyrosetta; print('PyRosetta OK')\" 2>&1 || echo 'PyRosetta: not available'")

    # ── Launch pipeline ────────────────────────────────────────────────────────
    print(f"\nLaunching pipeline (nohup, log → {REMOTE_LOG})...")
    exec_cmd(client, f"mkdir -p {REMOTE_RESULTS}/logs")

    run_cmd = (
        f"nohup python3 {REMOTE_SCRIPTS}/run_pipeline_direct.py "
        f"--data_dir {REMOTE_DATA} "
        f"--results_dir {REMOTE_RESULTS} "
        f"--parent_pdb /home/anees_22phd0670/8WT6.pdb "
        f"> {REMOTE_LOG} 2>&1 &"
    )
    print(f"  Command: {run_cmd}")
    exec_cmd(client, run_cmd)
    time.sleep(2)

    # Check it started
    rc, out = exec_cmd(client, "pgrep -a python3 | grep run_pipeline_direct")
    if out:
        print(f"  Pipeline PID: {out.split()[0]}")
        print(f"  Monitor: tail -f {REMOTE_LOG}")
    else:
        print("  WARNING: Pipeline process not found — check the log")

    exec_cmd(client, f"tail -20 {REMOTE_LOG} 2>/dev/null || echo 'Log not yet created'")

    sftp.close()
    client.close()
    print(f"\nDone. Pipeline running on VM.")
    print(f"Check progress: ssh {USERNAME}@{args.vm_host} 'tail -f {REMOTE_LOG}'")


if __name__ == "__main__":
    main()
