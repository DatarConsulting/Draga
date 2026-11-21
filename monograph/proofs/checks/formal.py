"""Build a pinned Lean project and inspect the kernel axiom report of every theorem.

The default project uses the standard library. The --real project also uses
the mathlib objects required for the real-valued Fibonacci constructions.
On a fresh checkout Lake retrieves the locked dependencies and those objects.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys

PROOFS = Path(__file__).resolve().parents[1]
CORE = {"propext", "Quot.sound", "Classical.choice"}


def cache_mathlib(root, source):
    imports = re.findall(r"^import (Mathlib\.[\w.]+)$", source, re.M)
    objects = root / ".lake/packages/mathlib/.lake/build/lib/lean"
    if all((objects / module.replace(".", "/")).with_suffix(".olean").is_file() for module in imports):
        return
    env = dict(os.environ, MATHLIB_NO_CACHE_ON_UPDATE="1")
    subprocess.run(["lake", "exe", "cache", "get", *imports], cwd=root,
                   env=env, stdout=sys.stderr, check=True)


def audit(real=False):
    directory, namespace = ("lean-real", "FibonacciMemory") if real else ("lean", "Draga")
    root = PROOFS / directory
    source = (root / (namespace + ".lean")).read_text()
    names = set(re.findall(r"^theorem (\w+)", source, re.M))
    declarations = re.findall(rf"^#print axioms {namespace}\.(\w+)$",
                              (root / "Audit.lean").read_text(), re.M)
    if not names or len(declarations) != len(names) or set(declarations) != names:
        raise ValueError("The axiom audit must cover every theorem exactly once")
    if re.search(r"^\s*(?:axiom|sorry|admit)\b", source, re.M):
        raise ValueError("Unverified declaration in Lean source")
    version = subprocess.check_output(["lean", "--version"], cwd=root, text=True).strip()
    if "version 4.19.0," not in version:
        raise ValueError("Use the pinned toolchain in " + directory + "/lean-toolchain")
    if real:
        cache_mathlib(root, source)
    subprocess.run(["lake", "build"], cwd=root, check=True, stdout=sys.stderr)
    raw = subprocess.check_output(["lake", "env", "lean", "Audit.lean"], cwd=root, text=True)
    dependencies = {}
    for name, tail in re.findall(rf"'{namespace}\.(\w+)' (.*)", raw):
        if tail == "does not depend on any axioms":
            axioms = []
        else:
            match = re.fullmatch(r"depends on axioms: \[(.*)\]", tail)
            if not match:
                raise ValueError("Unexpected Lean output " + tail)
            axioms = [s.strip() for s in match[1].split(",")]
        if name in dependencies or set(axioms) - CORE:
            raise ValueError(f"Duplicate report or unexpected axiom in {name}: {axioms}")
        dependencies[name] = axioms
    if set(dependencies) != names:
        raise ValueError("Incomplete kernel report")
    return {"toolchain": version, "theorems": len(names), "dependencies": dependencies}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--real", action="store_true", help="Audit the real-valued Fibonacci project")
    print(json.dumps(audit(parser.parse_args().real), indent=2))
