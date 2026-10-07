# Closing the Tsirelson Loophole in Relativistic Quantum Field Theory

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23202937.svg)](https://doi.org/10.5281/zenodo.23202937)

Code, data, figures and manuscript for

> **R. Chen**, *Closing the Tsirelson Loophole in Relativistic Quantum Field Theory: Split Inclusions, Hyperfinite
> Local Algebras, and the Effective Dimension of the Vacuum* (2026). DOI: [10.5281/zenodo.23202937](https://doi.org/10.5281/zenodo.23202937)

## Summary

MIP*=RE implies that commuting-operator correlations (C_qc) can lie at a finite distance from every limit of
finite-dimensional tensor-product correlations (C_qa). Such correlations would let a finite experiment detect an
actual infinity; this was the one loophole left by the certification budget of the companion paper
([doi:10.5281/zenodo.23202129](https://doi.org/10.5281/zenodo.23202129)).

- **Split inclusions.** With the split property (implied by Buchholz–Wichmann nuclearity, true for free fields), Bell
  correlations of local measurements in strictly spacelike separated regions are tensor-product correlations, in
  C_qs ⊆ C_qa.
- **Hyperfinite local algebras.** If Alice's local von Neumann algebra is injective (hyperfinite, as shown by
  Buchholz–D'Antoni–Fredenhagen under a scaling-limit assumption), every correlation with commuting partners lies in
  C_qa, even without a gap between the regions. Self-contained proof via semidiscreteness: Alice's measurements are
  compressed through matrix algebras by ucp maps.
- **de Sitter space.** The type II₁ observer algebra of Chandrasekaran–Longo–Penington–Witten is a corner of a crossed
  product by ℝ, hence injective under the same assumption.
- **Effective dimension.** An N-round Bell test certifies at most the smooth max-entropy of the local state at
  smoothing ε ~ 1/N. For the vacuum of a lattice free scalar field this follows an area law (≤ 9 bits for 128 sites at
  ε = 10⁻⁶, against 553 bits from naive counting) and grows polylogarithmically in 1/ε (2.3 → 13.1 bits from
  ε = 10⁻¹ to 10⁻¹²).

## Repository structure

| Path | Contents |
|---|---|
| `paper/` | LaTeX source and compiled PDF of the manuscript |
| `code/qft_effective_dimension.py` | Script producing every number and figure in the paper |
| `figures/` | Figure as vector PDF (used by the paper) and PNG |
| `data/` | Numerical data as CSV (metadata in `#` header lines) |
| `results/` | Console output of the script (the numbers quoted in the paper) |

## Reproducing the results

Requirements: Python ≥ 3.10 and the packages in `requirements.txt`
(tested with Python 3.12.4, NumPy 1.26.4, SciPy 1.13.1, Matplotlib 3.8.4). Runtime is a few seconds.

```bash
pip install -r requirements.txt
python code/qft_effective_dimension.py      # add --show to display the figure
```

Run from the repository root; the figure goes to `figures/` and data to `data/` (override with `QE_FIG_DIR`,
`QE_DATA_DIR`).

| Data file | Content | Paper |
|---|---|---|
| `block_effective_dimension.csv` | entanglement entropy and log2 d(ε) for blocks of the harmonic-chain vacuum | §4, Table 1, Fig. 1 |
| `dimension_vs_eps.csv` | d(ε) for m = 0.01, l = 64, ε = 10⁻¹ … 10⁻¹², with the product bound | §4, Fig. 1 |

To rebuild the paper (pdfLaTeX, two passes):

```bash
cd paper
pdflatex Chen_2026_QFT_Tsirelson_Loophole.tex
pdflatex Chen_2026_QFT_Tsirelson_Loophole.tex
```

## Citation

```bibtex
@misc{Chen2026QFTTsirelson,
  author = {Chen, Ruqing},
  title  = {Closing the Tsirelson Loophole in Relativistic Quantum Field Theory: Split Inclusions, Hyperfinite
            Local Algebras, and the Effective Dimension of the Vacuum},
  year   = {2026},
  doi    = {10.5281/zenodo.23202937},
  url    = {https://doi.org/10.5281/zenodo.23202937}
}
```

## License

- **Code** (`code/`): [MIT License](LICENSE)
- **Manuscript, figures, data and results** (`paper/`, `figures/`, `data/`, `results/`):
  [CC BY 4.0](LICENSE-CC-BY-4.0.md)

## Contact

Ruqing Chen — GUT Geoservice Inc., Montreal — ruqing@hotmail.com
