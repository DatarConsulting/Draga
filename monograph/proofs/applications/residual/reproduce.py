#!/usr/bin/env python3
"""Rebuild and replay the matched-residual experiment, without model weights."""
import argparse,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def run(args,log):
    print('Running',log,flush=True)
    with (ROOT/log).open('w') as f:subprocess.run(args,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,check=True)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repeats',type=int,default=3);ap.add_argument('--refit',action='store_true');a=ap.parse_args()
    if a.refit:run([sys.executable,'fit.py'],'fit.log')
    flags=['g++','-std=c++17','-O3','-march=native','-ffp-contract=off','-fno-fast-math','-I','src']
    run(flags+['src/experiment.cpp','-o','experiment'],'build.log')
    run(['./experiment','.','--test'],'results_unit.json')
    run(['./experiment','.',str(a.repeats)],'benchmark_progress.log')
    run(flags+['src/guard_benchmark.cpp','-o','guard_benchmark'],'build_guard.log')
    run(['./guard_benchmark','.'],'guard_stdout.log')
    run([sys.executable,'verify_independent.py'],'independent_stdout.log')
    run([sys.executable,'analyze.py'],'analysis_stdout.log')
    print('Done: Hankel_Exact_Residual_Report.html and results/summary.json')
if __name__=='__main__':main()
