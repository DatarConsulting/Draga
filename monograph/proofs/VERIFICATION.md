# Verification record

**Edition verified on 21 November 2026.** The reference is the 289-page book at repository commit `2a124be910a642d5cb16d85e52d95c2dbd75ce30`. The [statement index](INDEX.md) identifies the mathematical evidence for each result. The [machine-readable evidence](checks/evidence.json) records the figures below, including the individual Lean axiom dependencies.

| Check | Result | Interpretation |
|---|---|---|
| Chapter reconstruction | 3 of 3 original SHA-256 values matched | Expanding the proof includes and removing their file terminators returns exactly the pre-split source bytes. |
| Mathematical statements | 206 indexed, 197 explicit proof blocks attached to 195 statements | Four further statements have derivations in the surrounding section. Seven point to detailed results proved elsewhere in the book. |
| Numbering and citations | All 206 numbered PDF destinations, all 194 statement labels, all 590 distinct chapter labels and all 33 cited bibliography keys checked | This establishes source and cross-reference integrity. Bibliographic key resolution does not verify an external theorem. |
| Clean typesetting | All 289 pages matched in extracted text, geometry and 72-dpi rendering. The outlines matched. | Dates and binary PDF metadata are outside the comparison. |
| Lean 4.19.0 | 81 declared theorems built and audited, including 53 for the real-valued Fibonacci system | The kernel accepted the declarations under the hypotheses in the source. Their transitive logical dependencies use only standard Lean principles. |
| Executable Fibonacci codec | 1,968 exhaustive small-state transitions and 11,520 larger-threshold transitions. All 15,576 exact readouts matched across 2,088 traces. | The independent reference compared state identifiers and reduced rational midpoints. The uniform theorem covers every finite input history. [Evidence](checks/fibonacci-evidence.json) |
| Recovered software | 187 distinct historical source files, headers, tables or assembly listings across 45 source releases | All source bytes matched their recorded digests. All 151 distinct Python sources compiled. The 17 C and 2 C++ translation units passed syntax checks with their source trees. |
| Historical entry points | 85 of 85 completed with zero exit status | Runs include exact certificates, symbolic algebra, finite enumeration and numerical diagnostics. [Per-program ledger](checks/replay-results.json) |
| Independent rational fixtures | 12,297 substochastic tables, 6,762 stopped leaves and 123,702 individual updates on 5,694 maps | Integer or fraction arithmetic, with an equality case for the block budget and endpoint-separation checks for thresholds 1 through 16. |
| Finite transition graphs | 72 policies, 1,893 reachable states counted across policies, 3,786 exact labelled transitions | The packed native result and signed integer update agreed on every reachable edge. Induction extends this equality to all finite binary words on these machines. |
| Trajectory diagnostics | 52,637,040 prefix or stream steps across the 72 policies | Independent replay reproduced each archived maximum reference error within the documented binary64 tolerances. |
| Application pilots | 2,160 regenerated rational-pilot result rows, 10 checkpoint-bound certificates, 50 benchmark records and 1,048,576 bitwise word comparisons | Separate independent checkers audited rational bounds, checkpoint provenance or four serialized residual frames, as applicable. |
| Companion article | The notebook reran all 8 checkpoint instances. All 14 members of the regenerated arXiv source archive matched byte for byte. | Only the timing field differed in the raw inference record. |

## Fibonacci extension replay

A clean export restored the locked Lean dependencies and rebuilt both proof projects. The command-line suite and the proof notebook checked the Fibonacci results and compared all 289 pages to the reference. The [ordinal trace driver](lean-real/Trace.lean) produced the state and midpoint values recorded in the [worked example](lean-real/README.md).

The companion notebook regenerated the figures, PDFs and a byte-identical 14-member arXiv ZIP using the archived measurements. The earlier release records the independent eight-instance checkpoint replay. During this extension, automatic execution review blocked the optional inference run after ONNX Runtime attempted a telemetry connection to Microsoft. The new verification record counts the archived-measurement reproduction.

## Scope of the mathematical verification

The proof files preserve every printed argument. Lean checks specific quantified components of the operator obstruction, block-probability bound, quotient realization, invariant cover, active-map code and Fibonacci memory. For the latter it proves the exact persistent capacity, all-history correctness of the executable ordinal update, and the two regimes for intermediate tolerances. The [formalization table](README.md#what-the-formal-code-proves) states the assumptions and conclusion of each group. For example, the radix proof constructs both inverse maps and proves that an injective encoder for all reachable, individually queryable maps needs at least $(v+1)^k$ distinct codes.

The 81 Lean declarations cover the results and components listed in the formalization table. The inventory counts 206 numbered results in the book. Infinite-dimensional operator arguments, analytic limits, and the use of published Hardy-space extremals are documented in the written proofs. A numerical optimizer finishing successfully or a finite grid passing establishes evidence for that computation only. No independent peer review is recorded here.

## Reproduction details

The 85 entry points were run in isolated copies. Long programs that exceeded the initial 160-second diagnostic window were rerun with a longer allowance. The cone-oracle audit recomputed 119 reference values at 90 decimal digits. The largest observed error was approximately 0.39634 of its global bound. The arithmetic-filter verifier built its documented shared library. A probability pilot recovered its checksum-verified text fixture from the associated source release. The rank-two construction regenerated all 36 model files and matched their original decoded JSON values.

One historical stream driver resumes saved receipts when available. The separate [reachable-state verifier](checks/reachable.py) compared all 72 transition graphs and packet hashes from scratch. It reran the original trajectory schedules and recovered each of the 72 maximum observed errors. The real-valued reference tolerances are relative $10^{-10}$ and absolute $10^{-12}$. Integer coordinates, packet hashes and reachable-state counts were compared exactly.

Changing one proof byte, omitting a proof input, or supplying an unresolved theorem pointer caused the source verifier to reject the altered copy. The negative controls in the rational pilot detected the omission of a coupling term and of a quantization term. The checkpoint verifier rejected a deliberately damaged certificate. The full optional optimization-search and training programs are included as historical source. The reported checkpoint audit verifies the saved weights and recorded runs. It does not retrain the checkpoints.

The proof notebook was executed from a clean staged export on Linux. Google Colab is an additional entry point. Run commands, dependency versions, archive extraction and the mathematical role of each program are in the [proof guide](README.md).
