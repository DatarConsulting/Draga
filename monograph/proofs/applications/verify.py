"""Authenticate both pilot source trees and rerun their independent checks.

The rational recurrence pilot regenerates its synthetic trajectories and
per-sequence statistics with the released experiment before checking them.
The other pilots read the supplied fixed fixtures. All writes go to temporary
directories. Use --integrity-only to inspect the archival bytes without PyTorch.
"""
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent


def require(ok, desc):
    if not ok:
        raise ValueError(desc)


def restore(key, record, target):
    shutil.copytree(HERE / key, target)
    for name, digest in record["sources"].items():
        f = target / name
        require(sha256(f.read_bytes()).hexdigest() == digest, name)
        if f.suffix == ".py":
            compile(f.read_bytes(), str(f), "exec")
    fixture = HERE / (key + "-fixtures.zip")
    require(sha256(fixture.read_bytes()).hexdigest() == record["fixture_archive_sha256"], fixture)
    with ZipFile(fixture) as z:
        require(set(z.namelist()) == set(record["fixtures"]), "Fixture inventory")
        for name, digest in record["fixtures"].items():
            require(not Path(name).is_absolute() and ".." not in Path(name).parts, name)
            blob = z.read(name)
            require(sha256(blob).hexdigest() == digest, name)
            f = target / name
            f.parent.mkdir(parents=True, exist_ok=True)
            f.write_bytes(blob)


def run_member(root, name, args=()):
    entry = root / name
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    env.pop("PYTHONOPTIMIZE", None)
    subprocess.run([sys.executable, str(entry), *args], cwd=entry.parent, env=env, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--integrity-only", action="store_true")
    parser.add_argument("--pilot", choices=["all", "rational", "checkpoint", "residual"], default="all")
    args = parser.parse_args()
    records = json.loads((HERE / "manifest.json").read_text())
    with TemporaryDirectory(prefix="draga-pilots-") as temp:
        for key, record in records.items():
            path = Path(temp) / key
            restore(key, record, path)
            if args.integrity_only:
                continue
            if key == "residual" and args.pilot in ("all", "residual"):
                (path / "results").mkdir(exist_ok=True)
                run_member(path, "verify_independent.py")
            if key == "memory" and args.pilot in ("all", "rational"):
                base = "prior_evidence/HN_Memory_P1/"
                run_member(path, base + "experiment.py")
                run_member(path, base + "check_exact.py")
            if key == "memory" and args.pilot in ("all", "checkpoint"):
                base = "prior_evidence/certificate_pilot/code/"
                run_member(path, base + "verify_certificate.py", ["--require-complete"])
                run_member(path, base + "audit_checks.py", ["--require-complete"])
    print("Pilot checks completed", flush=True)


if __name__ == "__main__":
    main()
