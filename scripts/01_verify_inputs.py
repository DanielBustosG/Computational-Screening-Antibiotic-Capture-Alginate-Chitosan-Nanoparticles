"""01 - Verify raw inputs and report the analysed data snapshot.

Inputs : data/raw/md/*.xlsx, data/raw/docking/*.xlsx
Outputs: data/processed/input_manifest.csv

Checksums are recorded so that a later re-run can confirm it analysed the same
snapshot. A mismatch means the inputs changed, not that the analysis is wrong.
"""
import hashlib
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from algchi_screen.config import DATA_PROC, DATA_RAW


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    DATA_PROC.mkdir(parents=True, exist_ok=True)
    files = sorted(p for p in DATA_RAW.rglob("*") if p.is_file())
    if not files:
        raise SystemExit("No raw input files found under data/raw/.")
    manifest = pd.DataFrame([
        {"path": str(p.relative_to(DATA_RAW)),
         "size_bytes": p.stat().st_size,
         "sha256": sha256(p)}
        for p in files
    ])
    out = DATA_PROC / "input_manifest.csv"
    manifest.to_csv(out, index=False)
    print(f"files: {len(manifest)}")
    print(f"total size (MB): {manifest.size_bytes.sum() / 1e6:.1f}")
    print(f"manifest: {out.relative_to(DATA_PROC.parents[1])}")


if __name__ == "__main__":
    main()
