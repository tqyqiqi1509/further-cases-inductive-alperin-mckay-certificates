#!/usr/bin/env python3
"""Check the delivered file hashes; this does not execute the calculations."""
from pathlib import Path
import hashlib
root=Path(__file__).resolve().parent
count=0
for line in (root/"MANIFEST.sha256").read_text().splitlines():
    expected,name=line.split("  ",1)
    path=root/name
    if not path.is_file():raise SystemExit("MISSING: "+name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
        raise SystemExit("MISMATCH: "+name)
    count+=1
print(f"Verified {count} file hashes.")
