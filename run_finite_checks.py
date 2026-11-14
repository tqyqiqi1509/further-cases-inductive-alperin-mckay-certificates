#!/usr/bin/env python3
"""Reproduce the finite certificates; this does not verify the general proofs."""
from pathlib import Path
import argparse
import json
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
PYTHON_CHECKS = [
    "checks/check_a_hessian_exhaustive.py",
    "checks/check_root_models.py",
    "checks/check_d4_actual_root_images.py",
    "computations/f4_independent_integer.py",
    "computations/f4_small_field_source_round2.py",
    "computations/f4_small_field_parameters_round3.py",
]
GAP_CHECKS = [
    ("checks/o8_weyl_normalizers.g", "FRESH_O8_COMPLETE"),
    ("checks/o8_global_blocks.g", "FRESH_O8_GLOBAL_COMPLETE"),
]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="list checks without executing")
    parser.add_argument("--gap", help="GAP executable; otherwise use GAP on PATH")
    parser.add_argument("--gap-root", help="optional GAP root library, passed with -l")
    parser.add_argument("--python-only", action="store_true")
    args = parser.parse_args()
    if args.list:
        print("\n".join(PYTHON_CHECKS + [x[0] for x in GAP_CHECKS]))
        return 0

    out = ROOT / "reproduced_outputs"
    out.mkdir(exist_ok=True)
    jobs = [(p, [sys.executable, str(ROOT / p)], None) for p in PYTHON_CHECKS]
    if not args.python_only:
        gap = args.gap or shutil.which("gap")
        if not gap:
            parser.error("GAP is required, or choose --python-only for the Python checks")
        options = [gap, "-q", "-b"]
        if args.gap_root:
            options += ["-l", args.gap_root]
        jobs += [(p, options + [str(ROOT / p)], marker) for p, marker in GAP_CHECKS]

    report = []
    for name, command, marker in jobs:
        print(f"Running {name}", flush=True)
        result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
        stem = Path(name).stem
        (out / f"{stem}.stdout.txt").write_text(result.stdout)
        (out / f"{stem}.stderr.txt").write_text(result.stderr)
        errors = "Error," in result.stderr or "Syntax error" in result.stderr
        passed = result.returncode == 0 and not errors
        if marker:
            passed = passed and marker in result.stdout
        report.append({"script": name, "command": command,
                       "exit_code": result.returncode, "passed": passed})
        (out / "run_report.json").write_text(json.dumps(report, indent=2) + "\n")
        print("completed" if passed else "FAILED", flush=True)
        if not passed:
            return 1
    print("Finite checks completed. The general character and block proofs remain separate.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
