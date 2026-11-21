"""Run historical proof diagnostics in independent, disposable release trees.

Use --match to select an origin or script. A successful exit records that the
program finished with assertions enabled. The mathematical force of the run
depends on the propositions tested in that program.
"""
import argparse
from hashlib import sha256
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile
import time
from archive_tools import Catalogue

HERE = Path(__file__).resolve().parent


def replay(item, timeout):
    start = time.monotonic()
    with TemporaryDirectory(prefix="draga-replay-") as temp, Catalogue() as cat:
        root = Path(temp) / "release"
        cat.extract(item["origin"], root)
        entry = root / item["entry"]
        env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
                   MPLBACKEND="Agg", PYTHONDONTWRITEBYTECODE="1")
        env.pop("PYTHONOPTIMIZE", None)
        if entry.name == "symbolic_k3_response.py":
            (entry.parent / "research_v4_numeric").mkdir(exist_ok=True)
        if item["origin"] == "P4_Memoire_Lineaire_Preuve_Code" and entry.name == "verify.py":
            # The historical README builds this shared library before invoking its driver.
            subprocess.run(["cc", "-std=c11", "-O2", "-shared", "-fPIC", "arithmetic_filter.c",
                            "-o", "libarithmetic.so"], cwd=entry.parent, check=True)
        if entry.name == "extensions.py" and (entry.parent / "kernel.c").exists():
            subprocess.run(["cc","-std=c11","-O2","-shared","-fPIC","kernel.c","-o","kernel.so"],
                           cwd=entry.parent,check=True)
        if entry.name == "probability.py":
            # This text fixture was distributed with the associated memory-pilot release.
            name = "prior_evidence/certificate_pilot/data/pg15058_raw.txt"
            pilot = HERE.parent / "applications"
            expected = json.loads((pilot / "manifest.json").read_text())["memory"]["fixtures"][name]
            with ZipFile(pilot / "memory-fixtures.zip") as fixture:
                value = fixture.read(name)
            if sha256(value).hexdigest() != expected:
                raise ValueError("Text fixture digest mismatch")
            destination = entry.parent.parent / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(value)
        expected_models = {}
        if entry.name == "experiment.py" and item["origin"] == "P4_LaTeX_Preuve_Relecture_Complet":
            # The historical driver resumes cached candidates. Regenerate all 36 in this private copy.
            for cached in (entry.parent.parent / "models").glob("*.json"):
                expected_models[cached] = json.loads(cached.read_text())
                cached.unlink()
        process = subprocess.Popen([sys.executable, str(entry)], env=env, cwd=entry.parent,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        try:
            output, _ = process.communicate(timeout=timeout)
            status = "exit-zero" if process.returncode == 0 else "failed"
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            output, _ = process.communicate()
            status = "timeout"
        lines = output.decode(errors="replace").strip().splitlines()
        # A few historical diagnostics communicate failure via a JSON summary.
        negative = ('"all_passed": false', '"all_checks_passed": false', '"status": "fail"')
        if status == "exit-zero" and any(s in output.decode(errors="replace").lower() for s in negative):
            status = "reported-failure"
        extra = {}
        if status == "exit-zero" and expected_models:
            equal = sum(p.is_file() and json.loads(p.read_text()) == prior for p,prior in expected_models.items())
            extra = {"rebuilt_models_matching_original": equal, "models_expected": len(expected_models)}
            if equal != 36 or len(expected_models) != 36:
                status = "model-difference"
        return {**item, **extra, "status": status, "exit_code": process.returncode,
                "seconds": round(time.monotonic()-start, 2),
                "last_lines": [line[-600:] for line in lines[-5:]]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--match", default="")
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = json.loads((HERE / "entrypoints.json").read_text())
    selected = [r for r in rows if args.match in r["origin"]+"/"+r["entry"]]
    if not selected or min(args.jobs, args.timeout) < 1:
        parser.error("Select at least one entry point and positive limits.")
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(replay, item, args.timeout) for item in selected]
        results = []
        for job in futures:
            row = job.result()
            results.append(row)
            print(row["status"], row["origin"], Path(row["entry"]).name, flush=True)
    args.output.write_text(json.dumps(results, indent=2) + "\n")
    return int(any(r["status"] != "exit-zero" for r in results))


if __name__ == "__main__":
    sys.exit(main())
