# Historical program index

The paths below open the recovered source bytes. The final column identifies a replay entry point, a helper or an archived generated input. The catalogue records the SHA-256 and every displayed path has a matching original release and member name. The replay process reconstructs the original relative paths before running a program.

[Release catalogue](README.md) · [Replay ledger](../checks/replay-results.json) · [SHA-256 catalogue](catalogue.json)


## AAK Binary SpectralGap Contraction Manuscript v1 0 Sources

| Source | Purpose or interface | Use |
|---|---|---|
| [rational_matrices.py](code/AAK_Binary_SpectralGap_Contraction_Manuscript_v1_0_Sources/aak_manuscript_sources/rational_matrices.py) | Exact rational matrix helpers (Python standard library). | Supporting program |
| [verify_binary_certificate.py](code/AAK_Binary_SpectralGap_Contraction_Manuscript_v1_0_Sources/aak_manuscript_sources/verify_binary_certificate.py) | Exact finite certificates for the binary fixed-rank divergent-ratio family. | Replay |
| [verify_contraction_lift.py](code/AAK_Binary_SpectralGap_Contraction_Manuscript_v1_0_Sources/aak_manuscript_sources/verify_contraction_lift.py) | Exact interior-word certificates for the failed row-contractive lift. | Replay |

## AAK Fixed Alphabet Rank Tractability v1 0 Sources

| Source | Purpose or interface | Use |
|---|---|---|
| [verify_exact_identities.py](code/AAK_Fixed_Alphabet_Rank_Tractability_v1_0_Sources/AAK_Fixed_Alphabet_Rank_Tractability_v1_0/verify_exact_identities.py) | Finite rational checks for the AAK fixed-(d,k) proof note. | Replay |

## DRAGA AAK Certified Progress v0 8 0 Sources

| Source | Purpose or interface | Use |
|---|---|---|
| [explore_sections.py](code/DRAGA_AAK_Certified_Progress_v0_8_0_Sources/DRAGA_AAK_Certified_Progress_v0_8_0/explore_sections.py) | Uncertified section exploration. | Supporting program |
| [search_general.py](code/DRAGA_AAK_Certified_Progress_v0_8_0_Sources/DRAGA_AAK_Certified_Progress_v0_8_0/search_general.py) | Numerical exploration, never a certificate of global optimality. | Supporting program |
| [sos_compressed.py](code/DRAGA_AAK_Certified_Progress_v0_8_0_Sources/DRAGA_AAK_Certified_Progress_v0_8_0/sos_compressed.py) | Order-two moment SDP on a necessary compressed Hankel rank condition. | Supporting program |
| [verify_all.py](code/DRAGA_AAK_Certified_Progress_v0_8_0_Sources/DRAGA_AAK_Certified_Progress_v0_8_0/verify_all.py) | Replay the lower, upper and relaxation-barrier certificates exactly. | Replay |

## DRAGA AAK Compact Global Model v0 9 0 Sources

| Source | Purpose or interface | Use |
|---|---|---|
| [certify_relaxation.py](code/DRAGA_AAK_Compact_Global_Model_v0_9_0_Sources/certify_relaxation.py) | Rational reconstruction/checker for this compact moment SDP. | Replay |
| [compact_model.py](code/DRAGA_AAK_Compact_Global_Model_v0_9_0_Sources/compact_model.py) | Exact compact all-word PMI model; real two-state competitors, no commutation. | Supporting program |
| [certify_sos.py](code/DRAGA_AAK_Compact_Global_Model_v0_9_0_Sources/prior_v08/certify_sos.py) | Rational PSD Gram factors plus an exact coefficient-error budget. | Replay |
| [polynomial_model.py](code/DRAGA_AAK_Compact_Global_Model_v0_9_0_Sources/prior_v08/polynomial_model.py) | Exact polynomial data for a necessary four-dimensional compression. | Supporting program |
| [run_relaxation.py](code/DRAGA_AAK_Compact_Global_Model_v0_9_0_Sources/run_relaxation.py) | Symmetry-exact, order-two moment relaxation of the compact PMI model. | Supporting program |
| [summarize_results.py](code/DRAGA_AAK_Compact_Global_Model_v0_9_0_Sources/summarize_results.py) | Select only an exactly checked lower certificate and update display data. | Supporting program |
| [verify_all.py](code/DRAGA_AAK_Compact_Global_Model_v0_9_0_Sources/verify_all.py) | Replay the prerequisite and selected new certificate using exact arithmetic. | Replay |
| [verify_compact.py](code/DRAGA_AAK_Compact_Global_Model_v0_9_0_Sources/verify_compact.py) | Exact rational calibration and structure checks for the compact PMI model. | Replay |

## DRAGA AAK Constructive Controlled Gap v0 6 0 Sources

| Source | Purpose or interface | Use |
|---|---|---|
| [verify_constructive_gap.py](code/DRAGA_AAK_Constructive_Controlled_Gap_v0_6_0_Sources/aak_note_v06/verify_constructive_gap.py) | Exact recognizer and arithmetic replay for DRAGA AAK v0.6.0. | Replay |

## DRAGA AAK Fixed Target Certificates v0 7 0 Sources

| Source | Purpose or interface | Use |
|---|---|---|
| [explore_candidate.py](code/DRAGA_AAK_Fixed_Target_Certificates_v0_7_0_Sources/aak_note_v07/explore_candidate.py) | Numerical discovery only; use verify_fixed_target.py for certification. | Supporting program |
| [verify_fixed_target.py](code/DRAGA_AAK_Fixed_Target_Certificates_v0_7_0_Sources/aak_note_v07/verify_fixed_target.py) | Exact finite certificates for DRAGA v0.7.0 (not a general SOS solver). | Replay |

## DRAGA AAK Homogeneous Sharp Constant v0 4 0 Sources

| Source | Purpose or interface | Use |
|---|---|---|
| [verify_homogeneous.py](code/DRAGA_AAK_Homogeneous_Sharp_Constant_v0_4_0_Sources/verify_homogeneous.py) | Exact finite checks for DRAGA AAK v0.4.0 (Python 3 + SymPy). | Replay |

## DRAGA AAK Mixed Length Exact Optimum v0 5 0 Sources

| Source | Purpose or interface | Use |
|---|---|---|
| [verify_exact_optimum.py](code/DRAGA_AAK_Mixed_Length_Exact_Optimum_v0_5_0_Sources/verify_exact_optimum.py) | Exact replay and constructive recognition for DRAGA AAK v0.5.0. | Replay |

## DRAGA AAK Obstruction v0 3 0 Sources

| Source | Purpose or interface | Use |
|---|---|---|
| [verify_aak.py](code/DRAGA_AAK_Obstruction_v0_3_0_Sources/DRAGA_AAK_Obstruction_v0_3_0/verify_aak.py) | Exact arithmetic replay for DRAGA AAK v0.3.0. Python 3 + SymPy only. | Replay |

## DRAGA Commutative AAK Reproduction

| Source | Purpose or interface | Use |
|---|---|---|
| [verify_commutative_aak_gap.py](code/DRAGA_Commutative_AAK_Reproduction/verify_commutative_aak_gap.py) | Exact certificates for DRAGA_Commutative_AAK_Gap.pdf. | Replay |

## DRAGA Reviewed Complete

| Source | Purpose or interface | Use |
|---|---|---|
| [verify_review_changes.py](code/DRAGA_Reviewed_Complete/verify_review_changes.py) | Diagnostic checks for the revised DRAGA manuscripts. | Replay |

## Hankel Agent Memory Audit and Verification

| Source | Purpose or interface | Use |
|---|---|---|
| [verify_unequal_addresses.py](code/Hankel_Agent_Memory_Audit_and_Verification/verify_unequal_addresses.py) | Exact finite tests of the audit's rational stopping-address construction. | Replay |

## Hankel All Budgets Proof Bundle

| Source | Purpose or interface | Use |
|---|---|---|
| [hankel_boundary_checks.py](code/Hankel_All_Budgets_Proof_Bundle/hankel_boundary_checks.py) | Reproducible finite-Stein diagnostics of the binary boundary theorem. | Replay |

## Hankel Arxiv Source and Verification

| Source | Purpose or interface | Use |
|---|---|---|
| [verify_complex_hessian.py](code/Hankel_Arxiv_Source_and_Verification/anc/verify_complex_hessian.py) | Rebuild the three-state complex Hessian from the manuscript's local chart. | Replay |
| [verify_explicit_results.py](code/Hankel_Arxiv_Source_and_Verification/anc/verify_explicit_results.py) | Independent consistency checks for the explicit families in the manuscript. | Replay |

## Hankel Finite Optima Proof Bundle

| Source | Purpose or interface | Use |
|---|---|---|
| [f1_global_checks.py](code/Hankel_Finite_Optima_Proof_Bundle/f1_global_checks.py) | High-precision diagnostics for the exact rank-one binary Hankel optimum. | Replay |
| [g1_subroutine_checks.py](code/Hankel_Finite_Optima_Proof_Bundle/g1_subroutine_checks.py) | Reproduce small exact-arithmetic checks for the G1 construction. | Replay |
| [k2_first_correction_checks.py](code/Hankel_Finite_Optima_Proof_Bundle/k2_first_correction_checks.py) | F2 diagnostics, not proofs or global certificates. Requires numpy/scipy. | Replay |

## Hankel Passive Response Proof Bundle

| Source | Purpose or interface | Use |
|---|---|---|
| [b4_exact_response.py](code/Hankel_Passive_Response_Proof_Bundle/checks/b4_exact_response.py) | Exact response of the degree-three Brevig--Saksman extremizer. | Replay |
| [envelope_k3.py](code/Hankel_Passive_Response_Proof_Bundle/checks/envelope_k3.py) | Search the finite-dimensional three-state envelope numerically over the parameter grid. | Supporting program |
| [generalization_checks.py](code/Hankel_Passive_Response_Proof_Bundle/checks/generalization_checks.py) | Exact algebra checks and high-precision diagnostics for the v4 generalization. | Replay |
| [symbolic_k3_response.py](code/Hankel_Passive_Response_Proof_Bundle/checks/symbolic_k3_response.py) | Differentiate the passive objective in local coordinates and factor its Hessian characteristic polynomial. | Replay |

## Hankel Research Results Bundle

| Source | Purpose or interface | Use |
|---|---|---|
| [hankel_checks.py](code/Hankel_Research_Results_Bundle/hankel_checks.py) | Floating-point diagnostics for the word-Hankel structural-oracle manuscript. | Replay |

## Hankel Unified Manuscript Source and Evidence

| Source | Purpose or interface | Use |
|---|---|---|
| [verify_synthesis.py](code/Hankel_Unified_Manuscript_Source_and_Evidence/Hankel_Unified_Manuscript/verify_synthesis.py) | Independent checks of displayed formulas; standard library only. | Replay |

## Hankel v1 2 Math Audit

| Source | Purpose or interface | Use |
|---|---|---|
| [verify_patch_v1_2.py](code/Hankel_v1_2_Math_Audit/verify_patch_v1_2.py) | Focused reproducible mathematical checks for the v1.2 patch audit. | Replay |

## Hankel v1 3 Source and Verification

| Source | Purpose or interface | Use |
|---|---|---|
| [verify_rank_two_v1_3.py](code/Hankel_v1_3_Source_and_Verification/verify_rank_two_v1_3.py) | Exact algebra and finite checks for the version-1.3 additions. | Replay |

## Hankel v1 4 Source and Verification

| Source | Purpose or interface | Use |
|---|---|---|
| [discovery_audit.inc](code/Hankel_v1_4_Source_and_Verification/revision_v1_4/discovery_audit.inc) | Generated model constants or tables included by the native implementation. | Native input |
| [fixed_basis_sampling.inc](code/Hankel_v1_4_Source_and_Verification/revision_v1_4/fixed_basis_sampling.inc) | Generated model constants or tables included by the native implementation. | Native input |
| [integrate_v1_4.py](code/Hankel_v1_4_Source_and_Verification/revision_v1_4/integrate_v1_4.py) | Rebuild the standalone v1.4 source from its preserved v1.3 baseline. | Supporting program |
| [positive_discovery.inc](code/Hankel_v1_4_Source_and_Verification/revision_v1_4/positive_discovery.inc) | Generated model constants or tables included by the native implementation. | Native input |
| [verify_discovery_v1_4.py](code/Hankel_v1_4_Source_and_Verification/revision_v1_4/verify_discovery_v1_4.py) | Reproducible algebraic checks for v1.4. Tests are not general proofs. | Replay |

## Hankel v1 5 Source and Verification

| Source | Purpose or interface | Use |
|---|---|---|
| [audit.inc](code/Hankel_v1_5_Source_and_Verification/revision_v1_5/audit.inc) | Generated model constants or tables included by the native implementation. | Native input |
| [codec.inc](code/Hankel_v1_5_Source_and_Verification/revision_v1_5/codec.inc) | Generated model constants or tables included by the native implementation. | Native input |
| [discovery.inc](code/Hankel_v1_5_Source_and_Verification/revision_v1_5/discovery.inc) | Generated model constants or tables included by the native implementation. | Native input |
| [integrate_v1_5.py](code/Hankel_v1_5_Source_and_Verification/revision_v1_5/integrate_v1_5.py) | Rebuild standalone v1.5 from its preserved v1.4 baseline. | Supporting program |
| [memory.inc](code/Hankel_v1_5_Source_and_Verification/revision_v1_5/memory.inc) | Generated model constants or tables included by the native implementation. | Native input |
| [open.inc](code/Hankel_v1_5_Source_and_Verification/revision_v1_5/open.inc) | Generated model constants or tables included by the native implementation. | Native input |
| [quantitative.inc](code/Hankel_v1_5_Source_and_Verification/revision_v1_5/quantitative.inc) | Generated model constants or tables included by the native implementation. | Native input |
| [transfer.inc](code/Hankel_v1_5_Source_and_Verification/revision_v1_5/transfer.inc) | Generated model constants or tables included by the native implementation. | Native input |
| [verify_closures_v1_5.py](code/Hankel_v1_5_Source_and_Verification/revision_v1_5/verify_closures_v1_5.py) | Finite checks supporting v1.5; not formal proofs or an LLM benchmark. | Replay |

## Hankel v1 6 Source and Verification

| Source | Purpose or interface | Use |
|---|---|---|
| [initial_state.inc](code/Hankel_v1_6_Source_and_Verification/revision_v1_6/initial_state.inc) | Generated model constants or tables included by the native implementation. | Native input |
| [integrate_v1_6.py](code/Hankel_v1_6_Source_and_Verification/revision_v1_6/integrate_v1_6.py) | Build standalone v1.6 from the preserved v1.5 source and new proofs. | Supporting program |
| [resolvent.inc](code/Hankel_v1_6_Source_and_Verification/revision_v1_6/resolvent.inc) | Generated model constants or tables included by the native implementation. | Native input |
| [transcript_audit.inc](code/Hankel_v1_6_Source_and_Verification/revision_v1_6/transcript_audit.inc) | Generated model constants or tables included by the native implementation. | Native input |
| [verify_extension_v1_6.py](code/Hankel_v1_6_Source_and_Verification/revision_v1_6/verify_extension_v1_6.py) | Exact algebra and finite numerical checks, not formal proof certification. | Replay |

## Hankel v1 7 Source and Verification

| Source | Purpose or interface | Use |
|---|---|---|
| [compiler.py](code/Hankel_v1_7_Source_and_Verification/Hankel_v1_7_Source_and_Verification/implementation_snapshot/compiler.py) | Deterministic compiler for supplied symbolic memory contracts. | Supporting program |
| [build_local_runtime.py](code/Hankel_v1_7_Source_and_Verification/Hankel_v1_7_Source_and_Verification/revision_v1_7/build_local_runtime.py) | Optional reproducible workaround for the present unindexed TeX runtime. | Supporting program |
| [integrate_v1_7.py](code/Hankel_v1_7_Source_and_Verification/Hankel_v1_7_Source_and_Verification/revision_v1_7/integrate_v1_7.py) | Rebuild v1.7 from the preserved v1.6 baseline and reviewed fragments. | Supporting program |
| [verify_dhmc_v1_7.py](code/Hankel_v1_7_Source_and_Verification/Hankel_v1_7_Source_and_Verification/verify_dhmc_v1_7.py) | Exact finite checks supporting v1.7; not a proof of its general claims. | Replay |

## P0 P4 Audit Sans Tau Preuve Code

| Source | Purpose or interface | Use |
|---|---|---|
| [audit_controls.py](code/P0_P4_Audit_Sans_Tau_Preuve_Code/p4_cone_audit/code/audit_controls.py) | Negative controls, fallback certificate and adaptive precision stress. | Replay |
| [audit_extension.py](code/P0_P4_Audit_Sans_Tau_Preuve_Code/p4_cone_audit/code/audit_extension.py) | Executed audit: exact references are confined to this verifier. | Replay |
| [build_report_data.py](code/P0_P4_Audit_Sans_Tau_Preuve_Code/p4_cone_audit/code/build_report_data.py) | Build report tables from completed JSON, with an explicit preview mode. | Supporting program |
| [cases.py](code/P0_P4_Audit_Sans_Tau_Preuve_Code/p4_cone_audit/code/cases.py) | Supply the named parameter cases used by the companion verification programs. | Supporting program |
| [certificate.py](code/P0_P4_Audit_Sans_Tau_Preuve_Code/p4_cone_audit/code/certificate.py) | Rational sufficient certificates for positive binary 2x2 instruments. | Supporting program |
| [check_boundary.py](code/P0_P4_Audit_Sans_Tau_Preuve_Code/p4_cone_audit/code/check_boundary.py) | A deterministic rounding-boundary fixture; exact stored source coefficients. | Replay |
| [check_oracle_reference.py](code/P0_P4_Audit_Sans_Tau_Preuve_Code/p4_cone_audit/code/check_oracle_reference.py) | High-precision independent numerical sanity checks, NOT formal certificates. | Replay |
| [check_products.py](code/P0_P4_Audit_Sans_Tau_Preuve_Code/p4_cone_audit/code/check_products.py) | Independent exact projective products, and relative-error induction replay. | Replay |
| [check_root_leaf.py](code/P0_P4_Audit_Sans_Tau_Preuve_Code/p4_cone_audit/code/check_root_leaf.py) | Admissible narrow invariant interval: the root itself is a stopping leaf. | Replay |
| [codec.py](code/P0_P4_Audit_Sans_Tau_Preuve_Code/p4_cone_audit/code/codec.py) | Autonomous integer ID; no transition table and no history cache. | Supporting program |
| [dyadic.py](code/P0_P4_Audit_Sans_Tau_Preuve_Code/p4_cone_audit/code/dyadic.py) | Outward integer intervals. The operational power uses integer square roots. | Supporting program |
| [model.py](code/P0_P4_Audit_Sans_Tau_Preuve_Code/p4_cone_audit/code/model.py) | Parametric bounded interval kernel. No history/cache is retained here. | Supporting program |
| [oracle.py](code/P0_P4_Audit_Sans_Tau_Preuve_Code/p4_cone_audit/code/oracle.py) | Streaming transfer oracle with model-derived rational certificates. | Supporting program |
| [verify.py](code/P0_P4_Audit_Sans_Tau_Preuve_Code/p4_cone_audit/code/verify.py) | Independent rational belief-coordinate checks. Enumeration is test-only. | Replay |

## P0 P4 Comparaison Directe Code

| Source | Purpose or interface | Use |
|---|---|---|
| [audit_direct_grid.py](code/P0_P4_Comparaison_Directe_Code/p4_next_step/code/audit_direct_grid.py) | Independent exact-rational checks for the direct-grid baseline and selector. | Replay |
| [audit_hybrid_branch.py](code/P0_P4_Comparaison_Directe_Code/p4_next_step/code/audit_hybrid_branch.py) | Short interface smoke checks where the selector actually chooses the codec. | Replay |
| [direct_grid.py](code/P0_P4_Comparaison_Directe_Code/p4_next_step/code/direct_grid.py) | Certified finite-state baseline; state is one grid index, without a tree. | Supporting program |

## P4 Codage Uniforme Preuve Code

| Source | Purpose or interface | Use |
|---|---|---|
| [check_proof.py](code/P4_Codage_Uniforme_Preuve_Code/P4_Codage_Uniforme/code/check_proof.py) | Evaluate the finite inequalities and constants appearing in this release of the P4 proof. | Replay |
| [check_rolling.py](code/P4_Codage_Uniforme_Preuve_Code/P4_Codage_Uniforme/code/check_rolling.py) | Check the rolling state update on the enumerated paths. | Replay |
| [check_stream.py](code/P4_Codage_Uniforme_Preuve_Code/P4_Codage_Uniforme/code/check_stream.py) | Compare stream outputs with the exact reference on the selected sequences. | Replay |
| [driver.c](code/P4_Codage_Uniforme_Preuve_Code/P4_Codage_Uniforme/code/driver.c) | C translation unit. Exercises the native implementation. | Native C |
| [export_models.py](code/P4_Codage_Uniforme_Preuve_Code/P4_Codage_Uniforme/code/export_models.py) | Write the finite-model tables consumed by the native stream implementation. | Supporting program |
| [fold23.c](code/P4_Codage_Uniforme_Preuve_Code/P4_Codage_Uniforme/code/fold23.c) | C translation unit. Implements the fold23 routines of this release. | Native C |
| [uniform.py](code/P4_Codage_Uniforme_Preuve_Code/P4_Codage_Uniforme/code/uniform.py) | Uniform, exact integer construction. Immutable DP data is explicitly charged. | Supporting program |
| [uniform_filter.c](code/P4_Codage_Uniforme_Preuve_Code/P4_Codage_Uniforme/code/uniform_filter.c) | C translation unit. Implements the uniform filter routines of this release. | Native C |
| [verify.py](code/P4_Codage_Uniforme_Preuve_Code/P4_Codage_Uniforme/code/verify.py) | Run the release-specific verification driver on its included model or certificate instances. | Replay |

## P4 Compilation Robuste Preuve Code

| Source | Purpose or interface | Use |
|---|---|---|
| [check_boundary.py](code/P4_Compilation_Robuste_Preuve_Code/p4_robust/code/check_boundary.py) | Replay a rounding-boundary fixture found by exact parameter bisection. | Replay |
| [check_products.py](code/P4_Compilation_Robuste_Preuve_Code/p4_robust/code/check_products.py) | Long products against an exact rational oracle, used only in tests. | Replay |
| [check_proof.py](code/P4_Compilation_Robuste_Preuve_Code/p4_robust/code/check_proof.py) | Exact rational verification of all constants in the accompanying proof. | Replay |
| [check_stream_memory.py](code/P4_Compilation_Robuste_Preuve_Code/p4_robust/code/check_stream_memory.py) | Finite streams at precisions beyond old C limit, with Python memory audit. | Replay |
| [codec.py](code/P4_Compilation_Robuste_Preuve_Code/p4_robust/code/codec.py) | Autonomous arithmetic suffix code. Persistent history is the ID alone. | Supporting program |
| [dyadic.py](code/P4_Compilation_Robuste_Preuve_Code/p4_robust/code/dyadic.py) | Outward integer intervals; inherited series with explicit remainder bounds. | Supporting program |
| [model.py](code/P4_Compilation_Robuste_Preuve_Code/p4_robust/code/model.py) | Bounded-precision positive matrix products. All decisions are integer. | Supporting program |
| [oracle.py](code/P4_Compilation_Robuste_Preuve_Code/p4_robust/code/oracle.py) | Streaming positive transfer sums with bounded-precision dual trajectories. | Supporting program |
| [report_tables.py](code/P4_Compilation_Robuste_Preuve_Code/p4_robust/code/report_tables.py) | Convert recorded numerical run data into tables for the historical report. | Supporting program |
| [reprofile_memory.py](code/P4_Compilation_Robuste_Preuve_Code/p4_robust/code/reprofile_memory.py) | Reprofile existing certified configurations after a full-depth warm-up. | Replay |
| [verify.py](code/P4_Compilation_Robuste_Preuve_Code/p4_robust/code/verify.py) | Exhaustive small graphs, independent rational belief-coordinate oracle. | Replay |

## P4 Identifiant Direct Preuve Code

| Source | Purpose or interface | Use |
|---|---|---|
| [benchmark.c](code/P4_Identifiant_Direct_Preuve_Code/P4_Identifiant_Direct/code/benchmark.c) | C translation unit. Measures native filter operations. | Native C |
| [check_portable.py](code/P4_Identifiant_Direct_Preuve_Code/P4_Identifiant_Direct/code/check_portable.py) | Compare the portable identifier implementation with the reference on the prescribed sample. | Replay |
| [direct_filter.c](code/P4_Identifiant_Direct_Preuve_Code/P4_Identifiant_Direct/code/direct_filter.c) | C translation unit. Implements the direct filter routines of this release. | Native C |
| [portable.py](code/P4_Identifiant_Direct_Preuve_Code/P4_Identifiant_Direct/code/portable.py) | Exact arbitrary-precision direct suffix filter. | Supporting program |
| [sanitizer_driver.c](code/P4_Identifiant_Direct_Preuve_Code/P4_Identifiant_Direct/code/sanitizer_driver.c) | C translation unit. Exercises the native implementation under runtime instrumentation. | Native C |
| [verify_direct.py](code/P4_Identifiant_Direct_Preuve_Code/P4_Identifiant_Direct/code/verify_direct.py) | Exact exhaustive certificate; enumeration occurs only in this verifier. | Replay |

## P4 LaTeX Preuve Relecture Complet

| Source | Purpose or interface | Use |
|---|---|---|
| [build_manuscript.py](code/P4_LaTeX_Preuve_Relecture_Complet/P4_Manuscrit/build_manuscript.py) | Assemble the historical P4 manuscript and its supporting section sources. | Supporting program |
| [render_qa.py](code/P4_LaTeX_Preuve_Relecture_Complet/P4_Manuscrit/render_qa.py) | Render the historical book for page-level inspection. | Supporting program |
| [check_local_lemmas.py](code/P4_LaTeX_Preuve_Relecture_Complet/P4_Manuscrit/review/check_local_lemmas.py) | Narrow, independently written rational checks for the P4 referee review. | Replay |
| [belief.c](code/P4_LaTeX_Preuve_Relecture_Complet/P4_Manuscrit/supporting_code/breakthrough/code/belief.c) | C translation unit. Implements the belief routines of this release. | Native C |
| [build_report.py](code/P4_LaTeX_Preuve_Relecture_Complet/P4_Manuscrit/supporting_code/breakthrough/code/build_report.py) | Build the French execution report from recorded results, without rerunning tests. | Supporting program |
| [controls.py](code/P4_LaTeX_Preuve_Relecture_Complet/P4_Manuscrit/supporting_code/breakthrough/code/controls.py) | Boundary obstructions and diagnostic floating-point baselines. | Replay |
| [experiment.py](code/P4_LaTeX_Preuve_Relecture_Complet/P4_Manuscrit/supporting_code/breakthrough/code/experiment.py) | Constructive rank-two finite-bit pilot. No training, no claim of global optimum. | Replay |
| [extensions.py](code/P4_LaTeX_Preuve_Relecture_Complet/P4_Manuscrit/supporting_code/breakthrough/code/extensions.py) | Exact finite-state separation lower bounds and k=3,4 examples. | Replay |
| [final_audit.py](code/P4_LaTeX_Preuve_Relecture_Complet/P4_Manuscrit/supporting_code/breakthrough/code/final_audit.py) | Cross-check deployed packets, target rank, and the local-smoothing KL bound. | Replay |
| [kernel.c](code/P4_LaTeX_Preuve_Relecture_Complet/P4_Manuscrit/supporting_code/breakthrough/code/kernel.c) | C translation unit. Implements the kernel routines of this release. | Native C |
| [kl_realization.py](code/P4_LaTeX_Preuve_Relecture_Complet/P4_Manuscrit/supporting_code/breakthrough/code/kl_realization.py) | Same-state local smoothing, exact C output, rational TV and KL upper bounds. | Replay |
| [output_audit.py](code/P4_LaTeX_Preuve_Relecture_Complet/P4_Manuscrit/supporting_code/breakthrough/code/output_audit.py) | Exact C output check over complete reachable graphs, plus table-state baseline. | Replay |
| [probability.py](code/P4_LaTeX_Preuve_Relecture_Complet/P4_Manuscrit/supporting_code/breakthrough/code/probability.py) | Positive, normalized finite-belief machines. All acceptance tests use Fraction. | Replay |
| [probability_output_audit.py](code/P4_LaTeX_Preuve_Relecture_Complet/P4_Manuscrit/supporting_code/breakthrough/code/probability_output_audit.py) | Audit the probability outputs against the reference model and specified precision. | Replay |
| [replay.py](code/P4_LaTeX_Preuve_Relecture_Complet/P4_Manuscrit/supporting_code/breakthrough/code/replay.py) | Execute the actual packed C recurrence; independently replay using integers. | Replay |
| [verify.py](code/P4_LaTeX_Preuve_Relecture_Complet/P4_Manuscrit/supporting_code/breakthrough/code/verify.py) | Standalone rational certificate checker. Does not import the producer. | Replay |
| [verify_probability.py](code/P4_LaTeX_Preuve_Relecture_Complet/P4_Manuscrit/supporting_code/breakthrough/code/verify_probability.py) | Independent rational acceptance checker for the positive machines. | Replay |

## P4 Loi Memoire Sans Table Preuve Code

| Source | Purpose or interface | Use |
|---|---|---|
| [benchmark.c](code/P4_Loi_Memoire_Sans_Table_Preuve_Code/p4_pressure/code/benchmark.c) | C translation unit. Measures native filter operations. | Native C |
| [certify_pressure.py](code/P4_Loi_Memoire_Sans_Table_Preuve_Code/p4_pressure/code/certify_pressure.py) | Exact dyadic interval certificate: L_s h - h has a uniform strict sign. | Replay |
| [check_family.py](code/P4_Loi_Memoire_Sans_Table_Preuve_Code/p4_pressure/code/check_family.py) | Evaluate the defining family and the bounds on its derivative. | Replay |
| [compact_filter.c](code/P4_Loi_Memoire_Sans_Table_Preuve_Code/p4_pressure/code/compact_filter.c) | C translation unit. Implements the compact filter routines of this release. | Native C |
| [sanitizer_driver.c](code/P4_Loi_Memoire_Sans_Table_Preuve_Code/p4_pressure/code/sanitizer_driver.c) | C translation unit. Exercises the native implementation under runtime instrumentation. | Native C |
| [verify_compact.py](code/P4_Loi_Memoire_Sans_Table_Preuve_Code/p4_pressure/code/verify_compact.py) | Verify the compact filter construction and its finite reference cases. | Replay |

## P4 Memoire Lineaire Preuve Code

| Source | Purpose or interface | Use |
|---|---|---|
| [arithmetic_filter.c](code/P4_Memoire_Lineaire_Preuve_Code/p4_arithmetic/code/arithmetic_filter.c) | C translation unit. Implements the arithmetic filter routines of this release. | Native C |
| [certify_fast.py](code/P4_Memoire_Lineaire_Preuve_Code/p4_arithmetic/code/certify_fast.py) | Whole-domain certificate for the exact C fixed-point expression. | Replay |
| [check_proof.py](code/P4_Memoire_Lineaire_Preuve_Code/p4_arithmetic/code/check_proof.py) | Evaluate the finite inequalities and constants appearing in this release of the P4 proof. | Replay |
| [check_root.py](code/P4_Memoire_Lineaire_Preuve_Code/p4_arithmetic/code/check_root.py) | Replay of the executed root-refinement test (roughly 40 s locally). | Replay |
| [check_universal.py](code/P4_Memoire_Lineaire_Preuve_Code/p4_arithmetic/code/check_universal.py) | Test the universal stream representation against the release-specific reference. | Replay |
| [driver.c](code/P4_Memoire_Lineaire_Preuve_Code/p4_arithmetic/code/driver.c) | C translation unit. Exercises the native implementation. | Native C |
| [fast.py](code/P4_Memoire_Lineaire_Preuve_Code/p4_arithmetic/code/fast.py) | Integer specification of the finite C probability oracle; no floating point. | Supporting program |
| [universal.py](code/P4_Memoire_Lineaire_Preuve_Code/p4_arithmetic/code/universal.py) | Certified arbitrary-p arithmetic code. Integers only; no count table. | Supporting program |
| [verify.py](code/P4_Memoire_Lineaire_Preuve_Code/p4_arithmetic/code/verify.py) | Exhaustive comparison with independent rational belief-coordinate oracle. | Replay |
| [oracle.py](code/P4_Memoire_Lineaire_Preuve_Code/p4_arithmetic/reference/oracle.py) | Implement the Bayesian update and node generation used by the precision oracle. | Supporting program |
| [arithmetic_filter.asm](code/P4_Memoire_Lineaire_Preuve_Code/p4_arithmetic/results/arithmetic_filter.asm) | Historical disassembly of the arithmetic filter for inspecting generated machine code. | Compiler output |

## P4 Optimum Exact Preuve Code

| Source | Purpose or interface | Use |
|---|---|---|
| [optimal_memory.c](code/P4_Optimum_Exact_Preuve_Code/P4_Optimum/code/optimal_memory.c) | C translation unit. Implements the optimal memory routines of this release. | Native C |
| [robust_rank_four.py](code/P4_Optimum_Exact_Preuve_Code/P4_Optimum/code/robust_rank_four.py) | Exact checks for the perturbation corollary and explicit rank-four sources. | Replay |
| [verify_and_run.py](code/P4_Optimum_Exact_Preuve_Code/P4_Optimum/code/verify_and_run.py) | Exact tests of the rank-two memory theorem and full C transition graphs. | Replay |

## P4 Oracle Polynomial Preuve Code

| Source | Purpose or interface | Use |
|---|---|---|
| [audit.py](code/P4_Oracle_Polynomial_Preuve_Code/P4_Oracle_Polynomial/code/audit.py) | Exact references and enumeration occur only in this test/benchmark file. | Replay |
| [build_report.py](code/P4_Oracle_Polynomial_Preuve_Code/P4_Oracle_Polynomial/code/build_report.py) | Build the historical report from its stored results. | Supporting program |
| [spectral_codec.py](code/P4_Oracle_Polynomial_Preuve_Code/P4_Oracle_Polynomial/code/spectral_codec.py) | Experimental integration of the independently accepted polynomial oracle. | Supporting program |
| [spectral_oracle.py](code/P4_Oracle_Polynomial_Preuve_Code/P4_Oracle_Polynomial/code/spectral_oracle.py) | Polynomial transfer oracle: floating candidates, exact interval acceptance. | Supporting program |

## P4 Prediction Conditionnelle Preuve Code

| Source | Purpose or interface | Use |
|---|---|---|
| [check_analytic.py](code/P4_Prediction_Conditionnelle_Preuve_Code/p4_conditional/code/check_analytic.py) | Check the explicit identities used in the conditional-prediction argument. | Replay |
| [packed_predictor.c](code/P4_Prediction_Conditionnelle_Preuve_Code/p4_conditional/code/packed_predictor.c) | C translation unit. Implements the packed predictor routines of this release. | Native C |
| [sanitizer_driver.c](code/P4_Prediction_Conditionnelle_Preuve_Code/p4_conditional/code/sanitizer_driver.c) | C translation unit. Exercises the native implementation under runtime instrumentation. | Native C |
| [synthesize.py](code/P4_Prediction_Conditionnelle_Preuve_Code/p4_conditional/code/synthesize.py) | Synthesize certificates, WITHOUT treating experimental success as proof. | Supporting program |
| [verify.py](code/P4_Prediction_Conditionnelle_Preuve_Code/p4_conditional/code/verify.py) | Independent certificate checker: normalized matrix multiplication + Fractions. | Replay |

## P4 Reduction Surcout Preuve Code

| Source | Purpose or interface | Use |
|---|---|---|
| [audit_hybrid.py](code/P4_Reduction_Surcout_Preuve_Code/p4_reduced/code/audit_hybrid.py) | Evaluate the hybrid filter against the finite reference on the specified paths. | Replay |
| [audit_reduction.py](code/P4_Reduction_Surcout_Preuve_Code/p4_reduced/code/audit_reduction.py) | Exact rational verifier. Full-tree enumeration and memoization are TEST ONLY. | Replay |
| [build_report_data.py](code/P4_Reduction_Surcout_Preuve_Code/p4_reduced/code/build_report_data.py) | Collect recorded model and counterexample measurements for the historical report. | Supporting program |
| [compile_comparison.py](code/P4_Reduction_Surcout_Preuve_Code/p4_reduced/code/compile_comparison.py) | Compile the comparison programs used in the reduction-overhead experiment. | Supporting program |
| [reduced_codec.py](code/P4_Reduction_Surcout_Preuve_Code/p4_reduced/code/reduced_codec.py) | Integer prefix codec certified on internal nodes; no tree precomputation. | Supporting program |

## P4 Theoreme Parametrique Preuve Code

| Source | Purpose or interface | Use |
|---|---|---|
| [certificate.py](code/P4_Theoreme_Parametrique_Preuve_Code/p4_parametric/code/certificate.py) | Rational sufficient certificates for positive binary 2x2 instruments. | Supporting program |
| [check_certificates.py](code/P4_Theoreme_Parametrique_Preuve_Code/p4_parametric/code/check_certificates.py) | Sufficient-condition rejections and an admissibility boundary scan. | Replay |
| [check_neighborhood.py](code/P4_Theoreme_Parametrique_Preuve_Code/p4_parametric/code/check_neighborhood.py) | Reuse one unchanged certificate for two distinct true models in its inner box. | Replay |
| [check_oracle_reference.py](code/P4_Theoreme_Parametrique_Preuve_Code/p4_parametric/code/check_oracle_reference.py) | High-precision independent numerical sanity checks, NOT formal certificates. | Replay |
| [check_streams.py](code/P4_Theoreme_Parametrique_Preuve_Code/p4_parametric/code/check_streams.py) | Longer histories, model immutability and warm-state Python memory profiles. | Replay |
| [oracle.py](code/P4_Theoreme_Parametrique_Preuve_Code/p4_parametric/code/oracle.py) | Streaming transfer oracle with model-derived rational certificates. | Supporting program |
| [report_tables.py](code/P4_Theoreme_Parametrique_Preuve_Code/p4_parametric/code/report_tables.py) | Build report tables directly from completed results; never invent measurements. | Supporting program |
