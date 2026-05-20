#!/usr/bin/env python3
"""
Step 12b Boltz-1 Tier 2 runner — pen-stack/design:0.1.0 (boltz 2.2.1)
SAFE VERSION v3: Batched processing + correct FASTA header format.

ROOT CAUSE OF VM HANG:
  Boltz's preprocessor spawns min(32, n_cpu) threads regardless of input size.
  Feeding 188 sequences at once → 32 simultaneous feature-extraction threads
  → 20–40 GB RAM spike → kernel OOM → entire VM freezes (only 2 GB swap).

FIX (RAM OOM — VM freeze):
  - Process BATCH_SIZE=1 sequence per boltz predict call.
  - Between calls, subprocess exits, freeing all RAM before next sequence starts.
  - Docker --memory=40g provides a hard OOM kill ceiling (exit 137 not VM freeze).

FIX (VRAM OOM — "ran out of memory, skipping batch"):
  - With BATCH_SIZE > 1, boltz loads feature tensors for ALL sequences into VRAM
    simultaneously. For 10 long chimeric sequences (~400-700aa): > 16GB VRAM used.
  - Fix: BATCH_SIZE=1. Each boltz call has only 1 sequence's tensors in VRAM.
  - Cost: 188 individual subprocess calls, ~9 min each, total ~28h runtime.

FIX (Docker shm bus error):
  - num_workers=0 disables PyTorch DataLoader multiprocessing (no shm IPC).
  - Plus --shm-size=4g on docker run for belt-and-suspenders safety.

FIX (Critical — "Invalid record id"):
  Boltz 2.2.1 FASTA parser ONLY accepts headers in the form:
      >CHAIN_ID|entity_type   (e.g. >A|protein)
  ANY other format (including plain ">seq_id") raises:
      ValueError: Invalid record id: <name>
  Fix: write one FASTA file per sequence, named {seq_id}.fasta, containing
      >A|protein
      SEQUENCE
  Then pass the DIRECTORY of these files to boltz predict.
  Output is at: predictions/{seq_id}/A/confidence_A_model_0.json

FIX (MSA requirement — "Missing MSA's in input"):
  Boltz-1 requires Multiple Sequence Alignments (MSAs) for structure prediction.
  Without --use_msa_server, ALL sequences fail at preprocessing:
      RuntimeError: Missing MSA's in input and --use_msa_server flag not set.
  Fix: add --use_msa_server to boltz predict command.
  This uses the MMSeqs2/ColabFold public API (api.colabfold.com) to compute MSAs.
  Natural sequences (Strategy B): ColabFold returns rich MSAs → high-quality predictions.
  Novel chimeras (Strategy A/C/D): ColabFold finds distant homologs or returns
    sparse MSAs; boltz still produces predictions (lower but non-zero quality).
  Runtime note: MSA generation adds ~30-60s per sequence; total runtime ~34h.

Additional safeguards:
  - DONE.txt guard at startup (prevents restart loop after completion)
    with auto-invalidation of bug-era artifacts (gate6_pass=0, results>0)
  - Thermal monitor: SIGTERM boltz predict if GPU >= TEMP_ABORT_C
  - Checkpoint CSV every CHECKPOINT_EVERY results parsed
  - Heartbeat JSON every monitoring cycle
  - Boltz skips existing output (idempotent; safe to restart mid-batch)

Input:  /input/step12b_boltz1_tier2_input.fasta
Output: /output/
  boltz_predictions/              -- boltz writes here (batched subdirs)
    batch_0001/                   -- per-batch boltz out_dir
      boltz_results_b0001_input/  -- boltz creates this
        predictions/
          {seq_id}/A/confidence_A_model_0.json
          {seq_id}/A/{seq_id}_model_0.pdb
  _batch/                         -- per-batch INPUT dirs (temp per-seq FASTAs)
    b0001_input/{seq_id}.fasta
  step12b_results.csv
  step12b_results.parquet (if pandas available)
  heartbeat.json
  DONE.txt
"""

import sys, time, json, csv, subprocess, signal, threading, shutil
from pathlib import Path

sys.stdout.reconfigure(line_buffering=True)

# ── Paths ─────────────────────────────────────────────────────────────────────
FASTA_PATH     = Path("/input/step12b_boltz1_tier2_input.fasta")
OUTPUT_DIR     = Path("/output")
BOLTZ_OUT      = OUTPUT_DIR / "boltz_predictions"
BATCH_DIR      = OUTPUT_DIR / "_batch"          # temp per-batch FASTA
RESULTS_CSV    = OUTPUT_DIR / "step12b_results.csv"
HEARTBEAT_FILE = OUTPUT_DIR / "heartbeat.json"
DONE_TXT       = OUTPUT_DIR / "DONE.txt"

# ── Tuning ────────────────────────────────────────────────────────────────────
BATCH_SIZE        = 1     # sequences per boltz predict call
                          # MUST be 1: with 10 seqs, combined VRAM for all feature tensors
                          # exceeds 16 GB → "ran out of memory, skipping batch" for long seqs.
                          # 1 seq at a time: clean VRAM per call, ~9 min/seq, ~28h total.
CHECKPOINT_EVERY  = 20    # flush CSV every N completed sequences
TEMP_ABORT_C      = 92    # GPU °C: SIGTERM boltz, checkpoint, exit(1)
                          # RTX A4000 max safe ~93°C; 92 gives 1°C headroom.
COOL_BELOW_C      = 72    # wait between batches until GPU drops to this temp
COOL_POLL_S       = 10    # seconds between temperature polls during cooling
MONITOR_INTERVAL  = 15    # seconds between monitoring scans
NUM_WORKERS       = 0     # boltz dataloader workers: 0 = main-process only
                          # (avoids Docker shm bus error with --shm-size default 64MB)

# ── Gate 6 thresholds ─────────────────────────────────────────────────────────
PLDDT_THRESHOLD = 70.0
PTM_THRESHOLD   = 0.50

print(f"[Boltz-1 runner] Start — {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
print(f"[Boltz-1 runner] Batch size: {BATCH_SIZE} | Temp abort: {TEMP_ABORT_C}C", flush=True)

# ── DONE.txt guard: exit immediately if already finished ───────────────────────
if DONE_TXT.exists():
    try:
        done_data  = json.loads(DONE_TXT.read_text())
        r_parsed   = done_data.get("results_parsed", 0)
        g_pass     = done_data.get("gate6_pass", -1)
        total_seqs_done = done_data.get("total_sequences", 0)
        is_invalid = False

        # Bug-era artifact #1: completed but all boltz calls failed (no real PDBs)
        #   results_parsed > 0 but gate6_pass = 0
        if r_parsed > 0 and g_pass == 0:
            is_invalid = True
            reason = f"results_parsed={r_parsed} but gate6_pass=0 (all boltz calls failed)"

        # Bug-era artifact #2: thermal abort / early exit wrote DONE.txt before any results
        #   results_parsed = 0 but total_sequences > 0
        elif r_parsed == 0 and total_seqs_done > 0:
            is_invalid = True
            reason = f"results_parsed=0 (run aborted before any sequence completed)"

        if is_invalid:
            print(f"[Boltz-1 runner] DONE.txt is a bug-era artifact ({reason}). "
                  f"Removing and continuing.", flush=True)
            DONE_TXT.unlink()
        else:
            print("[Boltz-1 runner] DONE.txt exists — job already complete. Exiting.",
                  flush=True)
            print(DONE_TXT.read_text(), flush=True)
            sys.exit(0)
    except Exception:
        # Corrupt or unreadable DONE.txt — treat as complete to avoid loops
        print("[Boltz-1 runner] DONE.txt exists (unreadable) — "
              "treating as complete. Exiting.", flush=True)
        sys.exit(0)

# ── Verify CLI ────────────────────────────────────────────────────────────────
r = subprocess.run(["boltz", "predict", "--help"], capture_output=True)
if r.returncode not in (0, 1):
    print("[Boltz-1 runner] FATAL: boltz CLI not found", flush=True)
    sys.exit(1)
ver = subprocess.run(["python3", "-c", "import boltz; print(boltz.__version__)"],
                     capture_output=True, text=True).stdout.strip()
print(f"[Boltz-1 runner] boltz {ver} ready", flush=True)

import torch
print(f"[Boltz-1 runner] torch {torch.__version__}, CUDA={torch.cuda.is_available()}", flush=True)


# ── Helpers ───────────────────────────────────────────────────────────────────
def parse_fasta(path):
    seqs = {}
    hdr, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if hdr:
                    seqs[hdr[1:].split()[0]] = ("".join(seq), hdr)
                hdr, seq = line, []
            else:
                seq.append(line)
    if hdr:
        seqs[hdr[1:].split()[0]] = ("".join(seq), hdr)
    return seqs


def get_gpu_temp():
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=temperature.gpu", "--format=csv,noheader,nounits"],
            timeout=5).decode().strip()
        return int(out.split("\n")[0].strip())
    except Exception:
        return -1


def get_strategy(seq_id):
    sid = seq_id.upper()
    if sid.startswith("A_"): return "A"
    if "PROTMPNN" in sid or (sid.startswith("D") and "_IS621_" in sid): return "D"
    if "DEIMM" in sid or "IS621_DEIMM" in sid: return "C"
    return "B"


def find_confidence_json(seq_id):
    """Search boltz_predictions for this sequence's confidence JSON.

    Boltz 2.2.1 output layout (one FASTA per sequence named {seq_id}.fasta):
        boltz_predictions/batch_XXXX/boltz_results_bXXXX_input/
            predictions/{seq_id}/A/confidence_A_model_0.json

    The key invariant: the JSON is always named confidence_A_model_0.json
    and lives under a directory whose parent's name is the seq_id.
    i.e.  p.name == "confidence_A_model_0.json"
          p.parent.name == "A"
          p.parent.parent.name == seq_id
    """
    for p in BOLTZ_OUT.rglob("confidence_A_model_0.json"):
        if p.parent.parent.name == seq_id:
            return p
    # Legacy fallback (old single-fasta format — keeps resume safe)
    for p in BOLTZ_OUT.rglob(f"confidence_{seq_id}*.json"):
        return p
    return None


def find_pdb(seq_id):
    """Search for prediction PDB.

    Boltz 2.2.1 layout: predictions/{seq_id}/A/{seq_id}_model_0.pdb
    The PDB file's parent.parent.name == seq_id is the reliable check.
    """
    for p in BOLTZ_OUT.rglob("*.pdb"):
        if p.parent.parent.name == seq_id:
            return p
    return None


def is_done(seq_id):
    return find_confidence_json(seq_id) is not None or find_pdb(seq_id) is not None


def parse_scores(seq_id):
    """Return (mean_plddt 0-100, ptm 0-1)."""
    j = find_confidence_json(seq_id)
    if j:
        try:
            d = json.loads(j.read_text())
            raw = d.get("plddt", d.get("complex_plddt", 0.0))
            if isinstance(raw, list):
                plddt = float(sum(raw) / len(raw)) if raw else 0.0
            else:
                plddt = float(raw)
            if 0.0 <= plddt <= 1.0:   # scale 0-1 → 0-100
                plddt *= 100.0
            cptm = d.get("chains_ptm", None)
            if isinstance(cptm, list):   cptm = cptm[0] if cptm else 0.0
            elif isinstance(cptm, dict): cptm = next(iter(cptm.values()), 0.0)
            ptm = float(d.get("ptm", d.get("complex_ptm", cptm or 0.0)))
            return plddt, ptm
        except Exception as e:
            print(f"[Boltz-1 runner] JSON parse error {seq_id}: {e}", flush=True)
    # Fallback: B-factors
    pdb = find_pdb(seq_id)
    if pdb:
        bfacs = []
        with open(pdb) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    try: bfacs.append(float(line[60:66]))
                    except ValueError: pass
        if bfacs:
            return sum(bfacs) / len(bfacs), 0.0
    return 0.0, 0.0


def flush_csv(rows):
    with open(RESULTS_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader(); w.writerows(rows)


def write_heartbeat(batch_n, n_batches, done_seqs, total_seqs, n_pass, temp):
    HEARTBEAT_FILE.write_text(json.dumps({
        "ts":         time.strftime("%Y-%m-%d %H:%M:%S"),
        "batch":      batch_n,
        "total_batches": n_batches,
        "processed":  done_seqs,
        "total":      total_seqs,
        "gate6_pass": n_pass,
        "gpu_temp_c": temp,
    }))


# ── Load all sequences ─────────────────────────────────────────────────────────
sequences = parse_fasta(FASTA_PATH)
total_seqs = len(sequences)
print(f"[Boltz-1 runner] {total_seqs} sequences loaded", flush=True)

BOLTZ_OUT.mkdir(parents=True, exist_ok=True)
BATCH_DIR.mkdir(parents=True, exist_ok=True)

# ── Load existing CSV (resume) ─────────────────────────────────────────────────
CSV_FIELDS = ["id", "strategy", "seq_len", "mean_plddt", "ptm",
              "pass_plddt", "pass_ptm", "pass_gate6", "gpu_temp_c",
              "batch", "elapsed_s"]
csv_rows = []
done_csv_ids = set()
if RESULTS_CSV.exists():
    with open(RESULTS_CSV) as f:
        csv_rows = list(csv.DictReader(f))
        done_csv_ids = {r["id"] for r in csv_rows}
    print(f"[Boltz-1 runner] Loaded {len(csv_rows)} existing CSV rows", flush=True)

# sequences not yet in CSV
todo = [(sid, seq) for sid, seq in sequences.items() if sid not in done_csv_ids]
# but also skip if boltz already produced output (in case CSV flush was interrupted)
todo = [(sid, seq) for sid, seq in todo if not is_done(sid)]
print(f"[Boltz-1 runner] {len(todo)} sequences to process "
      f"({len(sequences) - len(todo)} already done)", flush=True)

if not todo:
    print("[Boltz-1 runner] Nothing to do — writing DONE.txt", flush=True)
    # parse any missing from CSV
    for sid, (seq, _) in sequences.items():
        if sid in done_csv_ids: continue
        if is_done(sid):
            plddt, ptm = parse_scores(sid)
            pp, pt = plddt >= PLDDT_THRESHOLD, ptm >= PTM_THRESHOLD
            csv_rows.append({"id": sid, "strategy": get_strategy(sid),
                              "seq_len": len(seq), "mean_plddt": f"{plddt:.3f}",
                              "ptm": f"{ptm:.4f}", "pass_plddt": pp, "pass_ptm": pt,
                              "pass_gate6": pp and pt, "gpu_temp_c": -1,
                              "batch": "N/A", "elapsed_s": "N/A"})
    flush_csv(csv_rows)
    n_pass = sum(1 for r in csv_rows if str(r.get("pass_gate6", "")) == "True")
    summary = {"total_sequences": total_seqs, "gate6_pass": n_pass,
               "gate6_fail": total_seqs - n_pass,
               "gate6_pass_rate": round(n_pass / total_seqs, 4) if total_seqs else 0,
               "completed_at": time.strftime("%Y-%m-%d %H:%M:%S")}
    DONE_TXT.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)
    sys.exit(0)

# ── Batch processing ───────────────────────────────────────────────────────────
batches = [todo[i:i + BATCH_SIZE] for i in range(0, len(todo), BATCH_SIZE)]
n_batches = len(batches)
print(f"[Boltz-1 runner] {n_batches} batches of up to {BATCH_SIZE} sequences each", flush=True)

t_run = time.time()
new_since_ckpt = 0
abort_flag = False
current_proc = [None]   # so thermal monitor can reach it


def thermal_monitor():
    """Background thread: check GPU temp every MONITOR_INTERVAL seconds."""
    while not abort_event.is_set():
        time.sleep(MONITOR_INTERVAL)
        temp = get_gpu_temp()
        if temp >= TEMP_ABORT_C and current_proc[0] is not None:
            print(f"[Boltz-1 runner] THERMAL ABORT: GPU {temp}C >= {TEMP_ABORT_C}C. "
                  f"Sending SIGTERM.", flush=True)
            flush_csv(csv_rows)
            current_proc[0].send_signal(signal.SIGTERM)
            abort_event.set()
            return


abort_event = threading.Event()
monitor_thread = threading.Thread(target=thermal_monitor, daemon=True)
monitor_thread.start()

for batch_idx, batch in enumerate(batches, 1):
    if abort_event.is_set():
        break

    batch_ids = [sid for sid, _ in batch]
    temp = get_gpu_temp()
    print(f"\n[Boltz-1 runner] Batch {batch_idx}/{n_batches}: "
          f"{len(batch)} sequences | GPU={temp}C", flush=True)

    # ── Write per-sequence FASTA files ──────────────────────────────────────
    # CRITICAL: Boltz 2.2.1 requires header format ">CHAIN_ID|entity_type"
    # A plain ">seq_id" header raises "ValueError: Invalid record id".
    # Fix: one file per sequence named {seq_id}.fasta, header ">A|protein".
    # Boltz identifies the sequence by the *filename stem* (seq_id), not header.
    batch_input_dir = BATCH_DIR / f"b{batch_idx:04d}_input"
    batch_input_dir.mkdir(parents=True, exist_ok=True)
    for sid, seq in batch:
        seq_fasta = batch_input_dir / f"{sid}.fasta"
        if not seq_fasta.exists():   # skip if already written (resume safety)
            seq_fasta.write_text(f">A|protein\n{seq}\n", encoding="utf-8")
    print(f"[Boltz-1 runner] Per-seq FASTAs written to {batch_input_dir}", flush=True)

    # ── Inter-batch GPU cooling ──────────────────────────────────────────────
    # After GPU inference, temp spikes. Wait until it drops before next call.
    if batch_idx > 1:
        cool_t = get_gpu_temp()
        if cool_t > COOL_BELOW_C:
            print(f"[Boltz-1 runner] Cooling: GPU {cool_t}C > {COOL_BELOW_C}C — waiting...",
                  flush=True)
            while True:
                time.sleep(COOL_POLL_S)
                cool_t = get_gpu_temp()
                if cool_t <= COOL_BELOW_C or abort_event.is_set():
                    break
                print(f"[Boltz-1 runner] Still cooling: {cool_t}C", flush=True)
            print(f"[Boltz-1 runner] Cooled to {cool_t}C — proceeding", flush=True)

    # Per-batch output subdirectory (avoids cross-batch file collisions)
    batch_out = BOLTZ_OUT / f"batch_{batch_idx:04d}"
    batch_out.mkdir(parents=True, exist_ok=True)

    # Pass the DIRECTORY of per-seq FASTAs to boltz predict
    cmd = [
        "boltz", "predict", str(batch_input_dir),
        "--out_dir",           str(batch_out),
        "--accelerator",       "gpu",
        "--output_format",     "pdb",
        "--recycling_steps",   "3",
        "--diffusion_samples", "1",
        "--model",             "boltz1",
        "--no_kernels",
        "--num_workers",       str(NUM_WORKERS),
        "--use_msa_server",    # required: boltz-1 needs MSAs; uses MMSeqs2/ColabFold API
    ]
    print(f"[Boltz-1 runner] Running: boltz predict {batch_input_dir.name} "
          f"--out_dir {batch_out.name} --model boltz1 --no_kernels", flush=True)

    t_batch = time.time()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, bufsize=1)
    current_proc[0] = proc

    for line in proc.stdout:
        line = line.rstrip()
        if line:
            print(f"[boltz] {line}", flush=True)

    proc.wait()
    current_proc[0] = None
    batch_elapsed = time.time() - t_batch
    print(f"[Boltz-1 runner] Batch {batch_idx} finished in {batch_elapsed:.0f}s "
          f"(exit {proc.returncode})", flush=True)

    if abort_event.is_set():
        break

    # Parse results for this batch
    temp = get_gpu_temp()
    for sid, seq in batch:
        if sid in done_csv_ids:
            continue
        plddt, ptm = parse_scores(sid) if is_done(sid) else (0.0, 0.0)
        pp = plddt >= PLDDT_THRESHOLD
        pt = ptm >= PTM_THRESHOLD
        row = {
            "id": sid, "strategy": get_strategy(sid), "seq_len": len(seq),
            "mean_plddt": f"{plddt:.3f}", "ptm": f"{ptm:.4f}",
            "pass_plddt": pp, "pass_ptm": pt, "pass_gate6": pp and pt,
            "gpu_temp_c": temp, "batch": batch_idx,
            "elapsed_s": f"{batch_elapsed / len(batch):.1f}",
        }
        csv_rows.append(row)
        done_csv_ids.add(sid)
        new_since_ckpt += 1

        if not is_done(sid):
            print(f"[Boltz-1 runner] WARNING: no output for {sid}", flush=True)

    total_done = len(done_csv_ids)
    n_pass = sum(1 for r in csv_rows if str(r.get("pass_gate6", "")) == "True")
    elapsed_h = (time.time() - t_run) / 3600
    rate = total_done / max(time.time() - t_run, 1)
    eta_h = ((total_seqs - total_done) / rate) / 3600 if rate > 0 else 0

    print(f"[Boltz-1 runner] Progress: {total_done}/{total_seqs} | "
          f"gate6_pass={n_pass} | elapsed={elapsed_h:.1f}h | ETA={eta_h:.1f}h", flush=True)

    write_heartbeat(batch_idx, n_batches, total_done, total_seqs, n_pass, temp)

    # Checkpoint CSV
    if new_since_ckpt >= CHECKPOINT_EVERY:
        flush_csv(csv_rows)
        print(f"[Boltz-1 runner] Checkpoint: {total_done} done, {n_pass} pass gate6", flush=True)
        new_since_ckpt = 0

    # Clean up batch input dir (save disk; prediction output is kept)
    try:
        shutil.rmtree(batch_input_dir)
    except Exception as e:
        print(f"[Boltz-1 runner] Could not clean batch input dir: {e}", flush=True)

# ── Stop monitor thread ────────────────────────────────────────────────────────
abort_event.set()
monitor_thread.join(timeout=30)

# ── Final save ─────────────────────────────────────────────────────────────────
flush_csv(csv_rows)
print(f"[Boltz-1 runner] CSV saved: {len(csv_rows)} rows", flush=True)

try:
    import pandas as pd
    df = pd.DataFrame(csv_rows)
    df.to_parquet(OUTPUT_DIR / "step12b_results.parquet", index=False)
    print(f"[Boltz-1 runner] Parquet saved", flush=True)
except Exception as e:
    print(f"[Boltz-1 runner] Parquet skipped: {e}", flush=True)

if abort_event.is_set():
    # Thermal abort — do not write DONE.txt; container will restart and resume.
    print("[Boltz-1 runner] Exiting with status 1 (thermal abort). "
          "Container will restart and resume.", flush=True)
    sys.exit(1)

n_total = len(csv_rows)
n_pass  = sum(1 for r in csv_rows if str(r.get("pass_gate6", "")) == "True")
elapsed_h = (time.time() - t_run) / 3600

summary = {
    "total_sequences": total_seqs,
    "results_parsed":  n_total,
    "gate6_pass":      n_pass,
    "gate6_fail":      n_total - n_pass,
    "gate6_pass_rate": round(n_pass / n_total, 4) if n_total else 0,
    "plddt_threshold": PLDDT_THRESHOLD,
    "ptm_threshold":   PTM_THRESHOLD,
    "batch_size":      BATCH_SIZE,
    "elapsed_h":       round(elapsed_h, 2),
    "completed_at":    time.strftime("%Y-%m-%d %H:%M:%S"),
}
DONE_TXT.write_text(json.dumps(summary, indent=2))
print(f"\n[Boltz-1 runner] COMPLETE — {n_pass}/{n_total} pass Gate 6 "
      f"in {elapsed_h:.1f}h", flush=True)
print(json.dumps(summary, indent=2), flush=True)
