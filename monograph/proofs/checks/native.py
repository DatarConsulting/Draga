"""Type-check every distinct historical C translation unit with its local headers.

The generated .inc files and model tables are restored at their original relative
paths. External linking and run-time behaviour are exercised separately by the
historical drivers that build and invoke the models.
"""
import argparse
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from archive_tools import Catalogue

ROOT = Path(__file__).resolve().parents[1]


def check_sources(compiler):
    rows = json.loads((ROOT / "archive/catalogue.json").read_text())
    result = []
    with Catalogue() as cat:
        for row in rows:
            if Path(row["member"]).suffix != ".c":
                continue
            with TemporaryDirectory(prefix="draga-c-") as temp:
                release = Path(temp) / "release"
                cat.extract(row["origin"], release)
                source = release / row["member"]
                p = subprocess.run([compiler, "-std=gnu11", "-fsyntax-only", str(source)],
                                   cwd=source.parent, capture_output=True, text=True)
                result.append(dict(source=row["path"], returncode=p.returncode,
                                   diagnostics=p.stderr.replace(str(release)+"/","")))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cc", default="cc")
    args = parser.parse_args()
    rows = check_sources(args.cc)
    print(json.dumps(rows, indent=2))
    raise SystemExit(int(any(r["returncode"] for r in rows)))
