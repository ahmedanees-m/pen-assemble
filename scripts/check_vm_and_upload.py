#!/usr/bin/env python3
"""
Check VM connectivity and upload pipeline data when accessible.
Run this script manually when you want to connect and launch the pipeline.

Usage:
    python check_vm_and_upload.py [--vm_host IP]
"""
import argparse
import io
import json
import sys
import time
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

try:
    import paramiko
except ImportError:
    print("ERROR: pip install paramiko")
    sys.exit(1)

KEY_FILE  = r"C:\Users\ANEES AHMED\.ssh\google_compute_engine"
USERNAME  = "anees_22phd0670"
LOCAL_BASE = Path(r"G:\My Drive\PEN-STACK\PAPER_4\pen-assemble")

REMOTE_DATA    = "/home/anees_22phd0670/pen_pipeline_data"
REMOTE_SCRIPTS = "/home/anees_22phd0670/pen_pipeline_scripts"
REMOTE_RESULTS = "/home/anees_22phd0670/pen_pipeline_results"
REMOTE_LOG     = f"{REMOTE_RESULTS}/logs/pipeline_run.log"

VM_IPS = ["10.110.5.102", "10.30.158.35", "10.74.0.51", "10.167.4.92"]


def try_connect(host: str, timeout: int = 8) -> tuple:
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    for key_type in [paramiko.RSAKey, paramiko.Ed25519Key, paramiko.ECDSAKey]:
        try:
            key = key_type.from_private_key_file(KEY_FILE)
            client.connect(host, username=USERNAME, pkey=key, timeout=timeout)
            return client, key_type.__name__
        except paramiko.AuthenticationException:
            continue
        except Exception:
            return None, "network_error"
    return None, "auth_failed"


def exec_cmd(client, cmd, print_output=True):
    _, stdout, stderr = client.exec_command(cmd)
    out = stdout.read().decode('utf-8', 'replace').strip()
    err = stderr.read().decode('utf-8', 'replace').strip()
    if print_output and out:
        print(f"  {out}")
    return out, err


def sftp_upload_dir(sftp, local_dir: Path, remote_dir: str):
    try: sftp.mkdir(remote_dir)
    except: pass
    for item in sorted(local_dir.iterdir()):
        remote_path = f"{remote_dir}/{item.name}"
        if item.is_dir():
            sftp_upload_dir(sftp, item, remote_path)
        else:
            size_kb = item.stat().st_size // 1024
            print(f"  Uploading {item.name} ({size_kb} KB) -> {remote_path}")
            sftp.put(str(item), remote_path)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--vm_host", default=None, help="VM IP (tries all known IPs if not set)")
    p.add_argument("--check_only", action="store_true", help="Only check connectivity, don't upload")
    p.add_argument("--status_only", action="store_true", help="Only check pipeline status on VM")
    args = p.parse_args()

    ips_to_try = [args.vm_host] if args.vm_host else VM_IPS

    # ── Find reachable VM ──────────────────────────────────────────────────────
    client = None
    connected_ip = None
    for ip in ips_to_try:
        print(f"Trying {ip}...", end=" ")
        c, info = try_connect(ip)
        if c:
            client = c
            connected_ip = ip
            print(f"CONNECTED ({info})")
            break
        else:
            print(f"failed ({info})")

    if not client:
        print("\nVM not reachable. Try again later.")
        print("Known IPs:", VM_IPS)
        sys.exit(1)

    # ── System info ────────────────────────────────────────────────────────────
    print("\nVM info:")
    exec_cmd(client, "hostname && uname -r")
    exec_cmd(client, "nvidia-smi --query-gpu=name,temperature.gpu,utilization.gpu --format=csv,noheader 2>/dev/null || echo 'no GPU info'")
    exec_cmd(client, "df -h /home 2>/dev/null | tail -1")

    if args.status_only:
        # Check pipeline status
        print("\nPipeline status:")
        exec_cmd(client, f"ls {REMOTE_RESULTS}/ 2>/dev/null && echo 'Results dir exists' || echo 'No results yet'")
        exec_cmd(client, f"ls {REMOTE_RESULTS}/step13_stability/ 2>/dev/null | head -5 || echo 'Step 13 not done'")
        exec_cmd(client, f"ls {REMOTE_RESULTS}/step15_mechclass/ 2>/dev/null | head -5 || echo 'Step 15 not done'")
        exec_cmd(client, f"ls {REMOTE_RESULTS}/step16_penscore/ 2>/dev/null | head -5 || echo 'Step 16 not done'")
        exec_cmd(client, f"tail -30 {REMOTE_LOG} 2>/dev/null || echo 'No log yet'")
        exec_cmd(client, "pgrep -a python3 | grep run_pipeline || echo 'Pipeline not running'")
        client.close()
        return

    if args.check_only:
        client.close()
        return

    # ── Check if pipeline already running ─────────────────────────────────────
    print("\nChecking if pipeline already running...")
    out, _ = exec_cmd(client, "pgrep -a python3 | grep run_pipeline_direct")
    if out:
        print(f"  Pipeline already running (PID: {out.split()[0]})")
        print(f"  Monitor: ssh {USERNAME}@{connected_ip} 'tail -f {REMOTE_LOG}'")
        client.close()
        return

    # ── Check if results already exist ────────────────────────────────────────
    out, _ = exec_cmd(client, f"test -f {REMOTE_RESULTS}/step16_penscore/pen_score_summary.json && echo DONE || echo NOT_DONE", print_output=False)
    if out.strip() == "DONE":
        print("\nPipeline already completed! Fetching results...")
        exec_cmd(client, f"cat {REMOTE_RESULTS}/step16_penscore/pen_score_summary.json")
        client.close()
        return

    # ── Upload data ────────────────────────────────────────────────────────────
    print(f"\nUploading pipeline data to VM ({connected_ip})...")
    sftp = client.open_sftp()

    # Check if data already uploaded
    out, _ = exec_cmd(client, f"ls {REMOTE_DATA}/designs/ 2>/dev/null | wc -l", print_output=False)
    if int(out.strip() or "0") >= 4:
        print("  Pipeline data already uploaded, skipping...")
    else:
        sftp_upload_dir(sftp, LOCAL_BASE / "pipeline_data", REMOTE_DATA)

    # Upload runner script
    script_local = LOCAL_BASE / "scripts/run_pipeline_direct.py"
    try: sftp.mkdir(REMOTE_SCRIPTS)
    except: pass
    sftp.put(str(script_local), f"{REMOTE_SCRIPTS}/run_pipeline_direct.py")
    print(f"  Uploaded run_pipeline_direct.py")
    sftp.close()

    # ── Check dependencies ─────────────────────────────────────────────────────
    print("\nChecking VM dependencies...")
    exec_cmd(client, "python3 -c \"import mech_class; print('mech_class:', mech_class.__version__)\" 2>&1")
    exec_cmd(client, "python3 -c \"import pen_score; print('pen_score:', pen_score.__version__)\" 2>&1")
    exec_cmd(client, "python3 -c \"import pyrosetta; print('PyRosetta OK')\" 2>&1 || echo 'PyRosetta: not available (Grantham fallback will be used)'")
    exec_cmd(client, f"test -f /home/{USERNAME}/8WT6.pdb && echo '8WT6.pdb: found' || echo '8WT6.pdb: NOT FOUND'")
    exec_cmd(client, f"ls /home/{USERNAME}/esm_tier1_output/pdbs/ | wc -l | xargs echo 'ESMFold PDBs:'")

    # ── Launch pipeline ────────────────────────────────────────────────────────
    print(f"\nLaunching pipeline...")
    exec_cmd(client, f"mkdir -p {REMOTE_RESULTS}/logs")

    nohup_cmd = (
        f"nohup python3 {REMOTE_SCRIPTS}/run_pipeline_direct.py "
        f"--data_dir {REMOTE_DATA} "
        f"--results_dir {REMOTE_RESULTS} "
        f"--parent_pdb /home/{USERNAME}/8WT6.pdb "
        f"> {REMOTE_LOG} 2>&1 &"
    )
    print(f"  {nohup_cmd}")
    exec_cmd(client, nohup_cmd)
    time.sleep(3)

    out, _ = exec_cmd(client, "pgrep -l python3 | grep run_pipeline_direct || echo 'NOT STARTED'")
    if "NOT STARTED" not in out:
        print(f"\nPipeline launched!")
        print(f"Monitor with: ssh {USERNAME}@{connected_ip} 'tail -f {REMOTE_LOG}'")
        exec_cmd(client, f"tail -20 {REMOTE_LOG} 2>/dev/null || echo 'Log not yet available'")
    else:
        print("\nWARNING: Pipeline may not have started. Check log manually.")
        exec_cmd(client, f"tail -20 {REMOTE_LOG} 2>/dev/null || echo 'No log'")

    # Save connection info for future use
    conn_info = {"ip": connected_ip, "username": USERNAME, "key_file": KEY_FILE,
                 "log": REMOTE_LOG, "results": REMOTE_RESULTS}
    (LOCAL_BASE / "pipeline_results_local_test/vm_connection.json").write_text(
        json.dumps(conn_info, indent=2)
    )

    client.close()
    print(f"\nDone. VM: {connected_ip}")


if __name__ == "__main__":
    main()
