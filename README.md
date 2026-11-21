# DRAGA

**Dimension Reduction and Guaranteed Approximation**  
Vashista Nobaub · Datar Consulting

[Read the monograph](monograph/monograph.pdf) · [Read the companion paper](paper/preprint.pdf) · [Read the mathematical supplement](paper/supplement.pdf) · [Reproduce the article in Colab](notebooks/DRAGA_Reproduce.ipynb) · [arXiv source archive](paper/arxiv-source.zip)

## Companion article

*When Nearly Rank-One Predictions Require Large Sequential Models* combines an operator-norm obstruction and an exact asymptotic state–log-loss tradeoff on one law of completed binary traces. The paper includes a CPU experiment on a pinned quantized Qwen2.5-0.5B-Instruct checkpoint. The analytic construction and the checkpoint observations have separate sections.

The self-contained English notebook runs without repository access after download. Its default mode regenerates the numerical tables, vector figures, article, standalone supplement and arXiv source archive. An optional cell downloads checksum-verified public weights and reruns all eight checkpoint instances. In Google Colab, open the notebook from this repository with your authorized account or upload the downloaded notebook, then use **Runtime → Run all**.

## Monograph

The book retains the three mathematical parts of the supplied 189-page November 2026 submission proof, with an editorial preface, navigation by topic and a companion-paper note. The book is typeset as a separate 289-page edition. All 590 chapter labels, 98 theorem environments and 197 proof environments from the original chapter bodies are retained. Broader context and original theorem numbering are preserved.

The [proof index](monograph/proofs/INDEX.md) links all 206 mathematical statements to their proofs or derivations. The [verification guide](monograph/proofs/README.md) includes the Lean formalization, the historical proof software and a separate [verification notebook](monograph/proofs/Verify_Monograph.ipynb). The [verification record](monograph/proofs/VERIFICATION.md) documents the scope of each check. The [Fibonacci proof and runnable codec](monograph/proofs/lean-real/README.md) establish the exact persistent-bit optimum for the two dyadic contractions, including the 90-bit instance.

With a standard TeX Live installation, run `latexmk -pdf monograph.tex` from `monograph/` or `latexmk -pdf main.tex` from `paper/`. All required figures are checked in. With Python, NumPy, SciPy and Matplotlib installed, run `python experiments/figures.py .` from the repository root to regenerate the figures and tables from the recorded inference data.

The arXiv source archive has a single document entry point, includes the proofs, and ships the inference records and Python programs as ancillary material. The official NeurIPS 2026 style is used in preprint mode. The style date identifies the template. Conference submission and editorial decisions are separate.
