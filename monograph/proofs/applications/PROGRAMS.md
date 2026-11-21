# Application program index

Every path below is checked against the SHA-256 recorded in the application manifest. The fixture packages support the independent checks described in the [pilot guide](README.md).

## Certified memory

| Program | Purpose or interface |
|---|---|
| [code/belief.c](memory/code/belief.c) | Native update for the discretized belief-state implementation. |
| [code/build_report.py](memory/code/build_report.py) | Build the French execution report from recorded results, without rerunning tests. |
| [code/controls.py](memory/code/controls.py) | Boundary obstructions and diagnostic floating-point baselines. |
| [code/experiment.py](memory/code/experiment.py) | Constructive rank-two finite-bit pilot. No training, no claim of global optimum. |
| [code/extensions.py](memory/code/extensions.py) | Exact finite-state separation lower bounds and k=3,4 examples. |
| [code/final_audit.py](memory/code/final_audit.py) | Cross-check deployed packets, target rank, and the local-smoothing KL bound. |
| [code/kernel.c](memory/code/kernel.c) | Packed fixed-point state transition used by the memory verifier. |
| [code/kl_realization.py](memory/code/kl_realization.py) | Same-state local smoothing, exact C output, rational TV and KL upper bounds. |
| [code/output_audit.py](memory/code/output_audit.py) | Exact C output check over complete reachable graphs, plus table-state baseline. |
| [code/probability.py](memory/code/probability.py) | Positive, normalized finite-belief machines. All acceptance tests use Fraction. |
| [code/probability_output_audit.py](memory/code/probability_output_audit.py) | Entry points `run`. |
| [code/replay.py](memory/code/replay.py) | Execute the actual packed C recurrence; independently replay using integers. |
| [code/verify.py](memory/code/verify.py) | Standalone rational certificate checker. Does not import the producer. |
| [code/verify_probability.py](memory/code/verify_probability.py) | Independent rational acceptance checker for the positive machines. |
| [prior_evidence/HN_Memory_P1/build_report.py](memory/prior_evidence/HN_Memory_P1/build_report.py) | Self-contained French report, with a compact theorem and measured scope. |
| [prior_evidence/HN_Memory_P1/check_exact.py](memory/prior_evidence/HN_Memory_P1/check_exact.py) | Independent rational replay, byte codec, and forward (all-controls) bounds. |
| [prior_evidence/HN_Memory_P1/codec.py](memory/prior_evidence/HN_Memory_P1/codec.py) | Fixed-width state codec. Only payload uint8 arrays are persistent memory. |
| [prior_evidence/HN_Memory_P1/experiment.py](memory/prior_evidence/HN_Memory_P1/experiment.py) | Approximate invariant flag, rational comparison certificate, held-out tests. |
| [prior_evidence/certificate_pilot/code/analyse.py](memory/prior_evidence/certificate_pilot/code/analyse.py) | Aggregate the frozen 50-run experiment; never train or select by test. |
| [prior_evidence/certificate_pilot/code/audit_checks.py](memory/prior_evidence/certificate_pilot/code/audit_checks.py) | Independent integrity audit, separate from model training and certification. |
| [prior_evidence/certificate_pilot/code/benchmark.py](memory/prior_evidence/certificate_pilot/code/benchmark.py) | Fresh-corpus pilot: scalar residual outputs and counting-stable multiletter heads. |
| [prior_evidence/certificate_pilot/code/certify.py](memory/prior_evidence/certificate_pilot/code/certify.py) | Exact dyadic stability and Hankel-error upper certificates for saved heads. |
| [prior_evidence/certificate_pilot/code/prepare_data.py](memory/prior_evidence/certificate_pilot/code/prepare_data.py) | Deterministic corpus preparation; no model output or test inspection. |
| [prior_evidence/certificate_pilot/code/replay_compression.py](memory/prior_evidence/certificate_pilot/code/replay_compression.py) | Prespecified rank8->4 replay. Numeric observations, not a roundoff proof. |
| [prior_evidence/certificate_pilot/code/run_all.py](memory/prior_evidence/certificate_pilot/code/run_all.py) | Entry points `run`. |
| [prior_evidence/certificate_pilot/code/verify_certificate.py](memory/prior_evidence/certificate_pilot/code/verify_certificate.py) | Independent verifier of serialized rational certificates. |

## Exact residual coding

| Program | Purpose or interface |
|---|---|
| [analyze.py](residual/analyze.py) | Entry points `main`. |
| [fit.py](residual/fit.py) | Calibration-only reconstruction of the archived finite Hankel selector. |
| [reference/fit_bases.py](residual/reference/fit_bases.py) | Calibration-only common value bases; weighted Stein proxy, not a Qwen certificate. |
| [reference/fit_finite_hankel.py](residual/reference/fit_finite_hankel.py) | Supplementary candidate selection using a complete finite scalar word-Hankel. |
| [reproduce.py](residual/reproduce.py) | Rebuild and replay the matched-residual experiment, without model weights. |
| [src/experiment.cpp](residual/src/experiment.cpp) | Driver for comparing the frozen predictors on the recorded states. |
| [src/gdn-local.hpp](residual/src/gdn-local.hpp) | C++ definitions used by the native prediction and encoding programs. |
| [src/gdn-lossless.hpp](residual/src/gdn-lossless.hpp) | C++ definitions used by the native prediction and encoding programs. |
| [src/gdn-predictive.hpp](residual/src/gdn-predictive.hpp) | C++ definitions used by the native prediction and encoding programs. |
| [src/guard_benchmark.cpp](residual/src/guard_benchmark.cpp) | Times the native guarded update on the selected model and input. |
| [src/independent_eval.c](residual/src/independent_eval.c) | Independent scalar fused multiply-add evaluator used by the packet verifier. |
| [verify_independent.py](residual/verify_independent.py) | Independent Python packet parser with a separate scalar FMA evaluator. |
