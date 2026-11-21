# Historical proof software

The [program index](PROGRAMS.md) covers 163 distinct source files, headers, tables and assembly listings. Their 325 recorded appearances span 43 source releases, including embedded archives. Identical SHA-256 values represent exactly the same bytes. Each distinct item also has a [readable working copy](code). The [catalogue](catalogue.json) records which original path supplied each working copy.

## Finding a computation

| Question | Where to look | What a replay establishes |
|---|---|---|
| How is the binary word-Hankel obstruction witnessed? | Binary fixed-rank and spectral-gap sources, `verify_binary_certificate.py` and `verify_contraction_lift.py` | Exact rational identities for the instances and lifted matrices tested. The quantifiers over competitors and the limiting argument are in §§1.8–1.12. |
| Which statements about boundary optima have symbolic checks? | Finite-optima and passive-response releases, especially `b4_exact_response.py`, `symbolic_k3_response.py` and `verify_complex_hessian.py` | Identities and factorizations of explicit matrices or scalar branches. A computed Hessian alone does not classify the full optimizing set. |
| Which computations concern stochastic realizations? | `verify_explicit_results.py` and the v1.2–v1.6 verifiers | Explicit moments, perturbation formulas, examples of failure and finite model diagnostics. Check the arithmetic and tolerances of the selected routine. |
| What happens when a history is replaced by a finite state? | `verify_dhmc_v1_7.py`, `verify_unequal_addresses.py` and the P4 releases | Reachable finite machines, decoding and update identities, finite certificates, error margins and selected long-stream measurements. |
| Where are the separate application pilots? | [Application evidence](../applications/README.md) | A rational recurrence, ten checkpoint-bound certificates and an independent bitwise residual parser have their own data and verification route. |

## Reconstructing a release

Run the following commands from `monograph/proofs`. The extraction target must be a new directory. List the complete names with the first command, then select a release or a script.

```bash
python checks/archive_tools.py
python checks/archive_tools.py --extract Hankel_Finite_Optima_Proof_Bundle /tmp/draga-finite-optima
python checks/replay.py --match verify_complex_hessian.py --output /tmp/draga-hessian.json
```

`sources.zip` contains the manifest and the source objects. The sixteen files under `support/` contain the support objects grouped by the first hexadecimal digit of their SHA-256. The extractor restores each selected path and checks its bytes before writing them. Documentation-only PDF renderings and precompiled host executables were omitted. Source, certificate data, reference tables and source-release digests were retained. The catalogue can reconstruct these selected members, while the original ZIP digest identifies the historical release from which they came.

Each replay runs with Python assertions enabled in a separate temporary release tree. The C arithmetic driver builds its shared library from source as specified in its original README. Timings come from the current machine. Some historical files are search or training drivers, formatting programs or reusable modules with no independent mathematical assertion. The program index makes their role visible.

