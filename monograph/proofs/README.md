# Proofs of the DRAGA monograph

This is the source of the proofs printed in the 289-page edition. Start with the [statement index](INDEX.md) to find an argument by number. Follow the [verification record](VERIFICATION.md) for the checks performed on this edition. The [verification notebook](Verify_Monograph.ipynb) runs from a repository checkout or an uploaded repository ZIP.

## Reading the source

The book has 98 theorems, 50 propositions, 44 lemmas and 14 corollaries. Within each chapter, `chapter.tex` holds the statements and the mathematical discussion between them. Each separately typeset proof is included from a file named after its statement. Two statements have two such files. Eleven other statements give the derivation in the surrounding section or refer to a more detailed result. Their index entries point to the relevant source lines or theorem numbers.

| Part | Mathematical object | Chapter source |
|---|---|---|
| I | Bounded word-Hankel operators, finite-rank competitors and passive realizations | [Multiletter approximation](operator-approximation/chapter.tex) |
| II | Proper laws on completed words, stochastic realizations and approximation error | [Probability](probability/chapter.tex) |
| III | Causal state, future observations, autonomous transitions and bit budgets | [Executable memory](executable-memory/chapter.tex) |

[inventory.json](inventory.json) records the owner, cross-references and SHA-256 of each explicit proof. The digest covers the bytes from the opening proof environment through its closing brace. Reinserting those bytes reproduces the three previously published chapter files. The final newline of the standalone fragment is the file terminator. The historical chapter fingerprints are fixed for this edition, so a future mathematical revision will be visible as a change to the reconstructed source.

## What the formal code proves

The [standard-library development](lean/Draga.lean) and the [real-valued Fibonacci development](lean-real/README.md) use Lean 4.19.0. The second uses a pinned mathlib revision. The kernel reports the logical dependencies of each of their 81 declared theorems. Several Lean declarations establish different implications of a single numbered result in the book.

| Book argument | Formal declarations | Domain and remaining assumptions |
|---|---|---|
| Th. 1.8.1, central-binomial obstruction | `vandermonde`, `coefficientMultiplicity`, `reflectBinomial`, `centralBinomial`, `rayleighCross` | Splitting binary words proves the central-binomial path identity for every degree. The cross-multiplied Rayleigh inequality is checked for arbitrary rational data. The annihilating polynomial and operator realization are in the written proof. |
| Chapter 1, word splitting and §1.12 | `shiftEntry`, `scalarSwap` | Words over an arbitrary alphabet. Scalar order symmetry follows from the stated multiplicative law and a commutative product. The bounded-operator factorization is in the written proof. |
| Th. 2.7.6, the block-mass inequality | `totalAdd`, `totalMono`, `blockBudget` | A finite list of states and blocks and nonnegative integer weights with row totals at most a common denominator. Division gives the rational case. The logarithmic loss limit also uses submultiplicativity, arithmetic–geometric mean and a two-state upper construction. |
| Prop. 3.10.2, exact quotient and packing | `replay`, `rightCongruence`, `distinguish`, `quotientRuns`, `quotientAnswers`, `minimalClasses`, `packingInjection` | All finite continuations, deterministic updates and arbitrary answer types. For approximate answers the packing premise requires that two histories admit no common answer after a distinguishing continuation. A metric gap greater than twice the tolerance gives this premise. |
| Prop. 3.10.3, invariant cover | `invariantCover`, `certifiedReadout` | A cell containing the initial state, closure under each transition and a query-uniform readout certificate on each cell. The conclusion holds at every time, including time zero. |
| Lemma 3.13.1, stopped leaves and exact covers | `fibonacciLeaves`, `endpoint_spacing`, `minimum_cover` | The two affine contractions give separated endpoints and matching upper and lower bounds for covers by translated intervals of the prescribed length. |
| Th. 3.13.2 and 3.13.4, exact executable memory | `zero_minimum`, `minimum_states`, `minimum_bits`, `ninety_bits` | Exact capacity for each startup in the leaf union and a packing bound for every startup in the unit interval. The tolerance is required at startup and after every finite history. |
| Prop. 3.13.3, executable ordinal | `ordinal_encode`, `ordinal_decode`, `ordinal_step_correct`, `midpoint_exact`, `ordinal_accuracy` | Rank and unrank are inverses. The concrete transition and rational decoder satisfy the uniform absolute-error bound. Bit-operation and temporary-space complexity are analyzed in the written proof. |
| Prop. 3.13.6, intermediate tolerances | `nondyadic_bounds` | Gives the exact optimum below three quarters of the dyadic interval length, and two-sided bounds from that threshold up to the next scale. |
| Prop. 3.15.9, the active map | `projectRun`, `activeMapFuture`, `mapsSeparated`, `rankBound`, `decodeEncode`, `encodeDecode`, `slotsLower`, `mapCapacity` | A deterministic factorization through the current map and any required auxiliary state. For a fixed key order, all radix words have a unique identifier. Full reachability and exact key queries supply the operational lower bound. The alternative sparse-code length is derived in the book. |

### Why the capacity lower bound needs a query contract

Let each of the $k$ keys take one of $v$ values or be absent. Fix the key order and put $b=v+1$. If the current digits are $d_0,\ldots,d_{k-1}$, the code is $\sum_{i=0}^{k-1}d_i b^i$. Repeated Euclidean division recovers the digits, including high-order zeros. Lean proves both inverse identities and the exact range $[0,b^k)$. A machine that must distinguish all these maps induces an injection into its available persistent states. Hence a width of $w$ fits all maps exactly when $2^w\ge b^k$. When a key query is unavailable or a map is unreachable, the injection premise must be established on the actual observable, reachable set. [Statement and proof](executable-memory/3.15.9.tex)

### Why the row budget appears in the KL proof

For a fixed block $B$, one state can maximize its probability, while a different state maximizes another block. Summing over blocks therefore costs as much as one row budget per state. The two rows $(1,0)$ and $(0,1)$ attain equality. The Lean inequality quantifies over every finite rational table after a common denominator has been cleared. For the stochastic matrices in the book, the same inequality over the reals follows by summing nonnegative entries, then applying substochasticity to each row. The rest of the lower bound requires counting distinct blocks and taking the limit in the stated order. [Complete proof](probability/2.7.6.tex)

### An external input to the Hardy bridge

The lower bound in Theorem 1.23.2 has its own proof. The matching upper construction uses a coefficient extremal. Brevig and Saksman determine the low-degree extremals used here. Brevig and Søvik supply the outer-extremal theorem for general degree. The formal and numerical programs in this directory do not establish those external theorems. [Written proof](operator-approximation/1.23.2.tex) · [Low-degree paper](https://arxiv.org/abs/1905.09547) · [General outer-extremal paper](https://arxiv.org/abs/2610.21788)

## Reproducing the checks

From this directory on Linux, use Python 3.12, a C compiler and the packages listed below. The notebook installs missing tools in a Colab runtime. The command-line checks operate without network access once dependencies are available.

```bash
python checks/verify.py
python checks/exact_examples.py
python checks/native.py
python checks/formal.py
python checks/formal.py --real
python checks/fibonacci_codec.py
OPENBLAS_NUM_THREADS=1 python checks/reachable.py
mkdir -p /tmp/draga-book-verification
(cd .. && latexmk -pdf -halt-on-error -outdir=/tmp/draga-book-verification monograph.tex)
python checks/verify.py --aux /tmp/draga-book-verification/monograph.aux \
  --pdf /tmp/draga-book-verification/monograph.pdf
python checks/replay.py --jobs 2 --timeout 3600 --output /tmp/draga-historical.json
```

To run the whole monograph suite with one command, use `python checks/verify_all.py /tmp/draga-proof-run`. The output directory must be new. Use `--pilot-python` to select a separate interpreter with PyTorch. The `--core` option runs the book, formal and exact checks plus pilot-file integrity.

The publication run used Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, SymPy 1.14.0, mpmath 1.3.0, PyMuPDF 1.26.6, GCC 14.2.0 and pdfTeX 1.40.26. Historic optimization-search drivers may additionally require CVXPY and a solver. Rechecking a supplied rational certificate is a different operation from searching for a new certificate.

The [archive catalogue](archive/README.md) explains the historical programs, where their inputs come from, and how to replay one entry point. The numerical diagnostics exercise examples and parameter grids. General analytic statements retain the hypotheses and deductions in their cited written proofs.

The [application checks](applications/README.md) also reproduce the independent verifiers for the numerical pilots used in the source dossier.


## A complete transition check on the finite pilot

For each of 36 model instances and both bit-allocation policies, the independent replay starts at the serialized initial state and explores both outgoing transitions of every reachable state. On each edge it compares the C status, packed bytes and signed coordinates with the arbitrary-precision reference update. An induction on word length then gives identical discrete execution for every finite binary word on these policies.

The continuous reference uses binary64. Its maximum observed error is checked against the original record on every binary prefix up to length 16 and six fixed sequences of length 100,000 for each policy. A matching floating-point maximum describes this specified experiment. The uniform error guarantee uses the mathematical certificate. [Reachability checker](checks/reachable.py) · [Scalar reference loop](checks/trace.c)
