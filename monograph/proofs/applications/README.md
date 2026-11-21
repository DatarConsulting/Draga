# Application evidence cited by the monograph's source dossier

Two empirical releases accompanied the mathematical sources. Their [38 source-file occurrences](PROGRAMS.md) are preserved in the original directory hierarchy. [Origin names, member hashes and fixture hashes](manifest.json) identify the precise versions. Several memory-pilot files also occur in the historical P4 collection.

| Pilot | What the independent code checks | Source |
|---|---|---|
| Rational recurrence | Exact replay of rational updates and codec bounds, including controls that expose an omitted coupling or a quantization term. The runner regenerates the synthetic trajectories from fixed seeds. | [Rational verifier](memory/prior_evidence/HN_Memory_P1/check_exact.py) |
| Checkpoint certificates | Reconstructs rational matrices from the serialized checkpoints, verifies exact certificate inequalities, checks outward rounding and exercises a deliberately damaged certificate. A separate integrity audit checks the fifty recorded runs. | [Certificate verifier](memory/prior_evidence/certificate_pilot/code/verify_certificate.py) · [Integrity audit](memory/prior_evidence/certificate_pilot/code/audit_checks.py) |
| Bitwise residual packets | A separate parser and scalar C evaluator reconstruct the serialized native state from the archived frames and compare the resulting 32-bit words exactly. This checks the supplied packets and arithmetic contract. | [Independent decoder](residual/verify_independent.py) · [C evaluator](residual/src/independent_eval.c) |

The archive stores every recovered program in these two releases, including analysis and fitting scripts. The checked-in fixtures cover the independent verifiers. Full predictor refitting uses the larger training inputs in the original residual experiment, identified by its digest in the manifest. All validation writes occur in a temporary directory.

```bash
python applications/verify.py --integrity-only
python applications/verify.py --pilot rational
python applications/verify.py --pilot residual
python applications/verify.py --pilot checkpoint
```

The checkpoint audit also requires PyTorch. Its recorded replay used PyTorch 2.6.0+cpu. The residual pilot requires GCC with C11 and scalar `fmaf` support. The underlying empirical comparison is about the specific frozen predictors and packets. Its primary rank-96 candidate used more bytes than its matched energy baseline.
