#!/usr/bin/env python3
"""
Step 12 ESMFold Tier 1 batch runner — pen-stack/plm:0.1.0 (esm-fold:v4) container
Hardened for unattended long-running operation:
  - Resumes from existing PDBs automatically
  - Dynamic chunk_size based on sequence length
  - CUDA OOM retry with smaller chunk, then skip
  - torch.cuda.empty_cache() after every sequence
  - Heartbeat file updated every sequence (watchdog-readable)
  - Checkpoint CSV flushed every CHECKPOINT_EVERY sequences

Gate 5 (P4 pre-registration): mean_pLDDT >= 70.0 AND pTM >= 0.50
"""

import os, sys, time, json, csv
from pathlib import Path

FASTA_PATH       = Path("/input/step12_esm_tier1_combined.fasta")
OUTPUT_DIR       = Path("/output")
PDB_DIR          = OUTPUT_DIR / "pdbs"
RESULTS_CSV      = OUTPUT_DIR / "step12_results.csv"
HEARTBEAT_FILE   = OUTPUT_DIR / "heartbeat.json"
DONE_TXT         = OUTPUT_DIR / "DONE.txt"
CHECKPOINT_EVERY = 10

PLDDT_THRESHOLD  = 70.0
PTM_THRESHOLD    = 0.50

sys.stdout.reconfigure(line_buffering=True)
print(f"[ESMFold runner] Start — {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)

import torch
print(f"[ESMFold runner] torch {torch.__version__}, CUDA={torch.cuda.is_available()}, "
      f"devices={torch.cuda.device_count()}", flush=True)
import esm
print(f"[ESMFold runner] fair-esm imported", flush=True)


def parse_fasta(path):
    seqs = []
    hdr, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if hdr is not None:
                    seqs.append((hdr[1:].split()[0], "".join(seq)))
                hdr, seq = line, []
            else:
                seq.append(line)
    if hdr is not None:
        seqs.append((hdr[1:].split()[0], "".join(seq)))
    return seqs


def chunk_size_for_length(seq_len):
    """Smaller chunks for longer sequences to control peak VRAM."""
    if seq_len <= 400:   return 128
    elif seq_len <= 600: return 64
    else:                return 32


sequences = parse_fasta(FASTA_PATH)
print(f"[ESMFold runner] Loaded {len(sequences)} sequences", flush=True)

PDB_DIR.mkdir(parents=True, exist_ok=True)
done_ids = {p.stem for p in PDB_DIR.glob("*.pdb")}
print(f"[ESMFold runner] Already done: {len(done_ids)} — will skip", flush=True)

# Load existing CSV rows for resume
csv_rows = []
if RESULTS_CSV.exists():
    with open(RESULTS_CSV) as f:
        reader = csv.DictReader(f)
        csv_rows = list(reader)
    print(f"[ESMFold runner] Loaded {len(csv_rows)} existing CSV rows", flush=True)

CSV_FIELDS = ["id","strategy","seq_len","mean_plddt","ptm",
              "pass_plddt","pass_ptm","pass_gate5","elapsed_s"]

def flush_csv(rows):
    with open(RESULTS_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

def write_heartbeat(processed, total, n_pass, last_id):
    HEARTBEAT_FILE.write_text(json.dumps({
        "ts": time.strftime("%Y-%m-%d %H:%M:%S"),
        "processed": processed,
        "total": total,
        "gate5_pass": n_pass,
        "last_id": last_id,
    }))

# Load model
print("[ESMFold runner] Loading ESMFold v1...", flush=True)
t0 = time.time()
model = esm.pretrained.esmfold_v1()
model = model.eval().cuda()
print(f"[ESMFold runner] Model loaded in {time.time()-t0:.1f}s", flush=True)
used  = torch.cuda.memory_allocated() / 1e9
total_vram = torch.cuda.get_device_properties(0).total_memory / 1e9
print(f"[ESMFold runner] VRAM: {used:.1f}/{total_vram:.1f} GB after model load", flush=True)


def infer_one(seq_id, sequence, chunk_size):
    model.set_chunk_size(chunk_size)
    with torch.no_grad():
        output = model.infer(sequence, num_recycles=3)
    mean_plddt = float(output["mean_plddt"].item())
    ptm        = float(output["ptm"].item())
    pdb_str    = model.output_to_pdb(output)[0]
    (PDB_DIR / f"{seq_id}.pdb").write_text(pdb_str)
    return mean_plddt, ptm


def get_strategy(seq_id):
    sid = seq_id.upper()
    if sid.startswith("A_"):  return "A"
    if sid.startswith("D0") or sid.startswith("D001") or \
       sid.startswith("D002") or "PROTMPNN" in sid or \
       (sid.startswith("D") and "_IS621_PROTMPNN" in sid): return "D"
    if "DEIMM" in sid or "IS621_DEIMM" in sid: return "C"
    return "B"   # UniProt accession


todo = [(sid, seq) for sid, seq in sequences if sid not in done_ids]
print(f"[ESMFold runner] {len(todo)} sequences to process", flush=True)

t_run = time.time()
processed = 0

for seq_id, seq in todo:
    t_seq = time.time()
    seq_len = len(seq)
    chunk   = chunk_size_for_length(seq_len)
    mean_plddt, ptm = 0.0, 0.0
    error_msg = ""

    try:
        mean_plddt, ptm = infer_one(seq_id, seq, chunk)
    except torch.cuda.OutOfMemoryError:
        torch.cuda.empty_cache()
        # Retry with halved chunk size
        smaller = max(16, chunk // 2)
        try:
            print(f"[ESMFold runner] OOM on {seq_id} (len={seq_len}, chunk={chunk}), "
                  f"retrying chunk={smaller}", flush=True)
            mean_plddt, ptm = infer_one(seq_id, seq, smaller)
        except Exception as e2:
            error_msg = f"OOM_SKIP:{e2}"
            print(f"[ESMFold runner] SKIP {seq_id} after OOM retry: {e2}", flush=True)
    except Exception as e:
        error_msg = str(e)[:80]
        print(f"[ESMFold runner] ERROR {seq_id}: {e}", flush=True)
    finally:
        torch.cuda.empty_cache()

    elapsed       = time.time() - t_seq
    pass_plddt    = mean_plddt >= PLDDT_THRESHOLD
    pass_ptm      = ptm >= PTM_THRESHOLD
    pass_gate5    = pass_plddt and pass_ptm
    strategy      = get_strategy(seq_id)
    processed    += 1

    row = {"id": seq_id, "strategy": strategy, "seq_len": seq_len,
           "mean_plddt": f"{mean_plddt:.3f}", "ptm": f"{ptm:.4f}",
           "pass_plddt": pass_plddt, "pass_ptm": pass_ptm,
           "pass_gate5": pass_gate5, "elapsed_s": f"{elapsed:.1f}"}
    csv_rows.append(row)

    elapsed_total = time.time() - t_run
    rate          = processed / elapsed_total
    remaining     = len(todo) - processed
    eta_h         = (remaining / rate) / 3600 if rate > 0 else 0
    n_pass        = sum(1 for r in csv_rows if str(r["pass_gate5"]) == "True")

    print(f"[ESMFold runner] [{processed}/{len(todo)}] {seq_id[:40]} len={seq_len} "
          f"chunk={chunk} | pLDDT={mean_plddt:.1f} pTM={ptm:.3f} "
          f"{'PASS' if pass_gate5 else 'fail'} | {elapsed:.0f}s | "
          f"ETA {eta_h:.1f}h | gate5={n_pass}", flush=True)

    write_heartbeat(processed, len(todo), n_pass, seq_id)

    if processed % CHECKPOINT_EVERY == 0:
        flush_csv(csv_rows)
        print(f"[ESMFold runner] Checkpoint @ {processed} ({n_pass} pass gate5)", flush=True)

# Final save
flush_csv(csv_rows)
try:
    import pandas as pd
    df = pd.DataFrame(csv_rows)
    df.to_parquet(OUTPUT_DIR / "step12_results.parquet", index=False)
    print(f"[ESMFold runner] Parquet saved: {len(df)} rows", flush=True)
except Exception as e:
    print(f"[ESMFold runner] Parquet failed (CSV is good): {e}", flush=True)

n_pass  = sum(1 for r in csv_rows if str(r["pass_gate5"]) == "True")
n_total = len(csv_rows)
elapsed_h = (time.time() - t_run) / 3600

summary = {
    "total_sequences": n_total, "gate5_pass": n_pass,
    "gate5_fail": n_total - n_pass,
    "gate5_pass_rate": round(n_pass / n_total, 4) if n_total else 0,
    "plddt_threshold": PLDDT_THRESHOLD, "ptm_threshold": PTM_THRESHOLD,
    "elapsed_h": round(elapsed_h, 2),
    "completed_at": time.strftime("%Y-%m-%d %H:%M:%S"),
}
DONE_TXT.write_text(json.dumps(summary, indent=2))
print(f"\n[ESMFold runner] COMPLETE — {n_pass}/{n_total} pass Gate 5 "
      f"({100*n_pass/n_total:.1f}%) in {elapsed_h:.1f}h", flush=True)
print(json.dumps(summary, indent=2), flush=True)
