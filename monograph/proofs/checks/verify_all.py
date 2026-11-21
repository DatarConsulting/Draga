"""Run the monograph release checks and write all fresh outputs to a new directory.

The default includes the formal build, TeX comparison, historical replays and
independent application pilots. --core omits the long historical and pilot
runs. Use --pilot-python when PyTorch is installed in a separate environment.
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys

PROOFS = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path,help="A new report directory")
    parser.add_argument("--core",action="store_true")
    parser.add_argument("--jobs",type=int,default=2)
    parser.add_argument("--timeout",type=int,default=3600)
    parser.add_argument("--pilot-python",default=sys.executable)
    args = parser.parse_args()
    if min(args.jobs,args.timeout)<1:
        parser.error("The worker count and time limit must be positive")
    output = args.output.resolve()
    output.mkdir(parents=True,exist_ok=False)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
               PYTHONDONTWRITEBYTECODE="1")
    env.pop("PYTHONOPTIMIZE",None)

    def run(name,command,where=PROOFS):
        with (output/name).open("w") as log:
            subprocess.run(command,cwd=where,env=env,stdout=log,check=True)
        print("Verified",name,flush=True)

    def check(script,*params):
        return [sys.executable,"checks/"+script+".py",*params]

    run("lean.json",check("formal"))
    run("fibonacci-lean.json", check("formal", "--real"))
    run("fibonacci-codec.json", check("fibonacci_codec"))
    run("exact.json",check("exact_examples"))
    run("native.json",check("native"))
    run("reachability.json",check("reachable"))
    run("pilot-integrity.log",[sys.executable,"applications/verify.py","--integrity-only"])
    tex = output/"tex"
    tex.mkdir()
    run("typesetting.log",["latexmk","-pdf","-halt-on-error","-interaction=nonstopmode",
                           "-outdir="+str(tex),"monograph.tex"],PROOFS.parent)
    run("corpus.json",check("verify","--aux",tex/"monograph.aux","--pdf",tex/"monograph.pdf"))
    if not args.core:
        run("historical.log",check("replay","--jobs",str(args.jobs),"--timeout",str(args.timeout),
                                   "--output",output/"historical.json"))
        for pilot in ("rational","residual","checkpoint"):
            interpreter = args.pilot_python if pilot=="checkpoint" else sys.executable
            run("pilot-"+pilot+".log",[interpreter,"applications/verify.py","--pilot",pilot])
    print(output)


if __name__ == "__main__":
    main()
