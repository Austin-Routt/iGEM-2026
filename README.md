# AKGut — Flux-Balance Analysis of α-Ketoglutarate Overproduction in *E. coli* Nissle 1917

**UH iGEM 2026 (Houston, TX) · Dry Lab**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Wiki: <https://2026.igem.wiki/houston-tx/model> · Mirror of <https://github.com/Austin-Routt/iGEM-2026>

---

## What this software does

This repository contains the computational strain-design work behind **AKGut**, our project to engineer the probiotic *Escherichia coli* Nissle 1917 (EcN) to secrete **α-ketoglutarate (α-KG)** — the metabolite that extends *C. elegans* lifespan, with maximal effect at 8 mM (Chin et al., 2014).

Using **flux balance analysis (FBA)** on the published EcN genome-scale model **iHM1533** (van 't Hof et al., 2022; 3,143 reactions / 2,261 metabolites / 1,533 genes), the code answers one question: *which genetic edits force α-KG out of the cell without stalling growth?* It compares three strains — wild-type, our computational design **LitOpt**, and the **wet-lab build** — on both glucose and glycerol, and regenerates every figure on our Model wiki page.

**Who it's for.** iGEM teams and metabolic engineers who want to (a) reproduce our published numbers exactly, (b) re-run the design under their own constraints (different carbon source, oxygen level, growth floor), or (c) reuse the strain-design scaffolding for a different target metabolite. Everything is plain Python with one command to reproduce.

## Main results this code reproduces

All values in mmol·gDW⁻¹·h⁻¹, at the operating point (50 % of the aerobic O₂ maximum, growth held ≥ 80 % of each strain's own maximum):

| Strain | Glucose | Glycerol |
|---|---|---|
| Wild-type | 4.88 | 2.01 |
| **LitOpt** (computational design) | **6.10** | **2.38** |
| Wet-lab build | 4.88 | 2.01 |

Plus four findings the figures make explicit:

1. **Oxygen optimum is non-monotonic.** α-KG peaks near 50 % O₂ (microaerobic), not at full aeration — O₂ baselines 21.97 (glucose) / 13.53 (glycerol).
2. **PPC and PYC are stoichiometrically equivalent** at the operating point (6.10 either way). Pyruvate carboxylase is preferred for a *kinetic* reason FBA cannot represent: glucose import via the PTS spends PEP, so overexpressing PEP carboxylase competes with sugar uptake. This independently matches the notebook's Stage 13 conclusion.
3. **No design is growth-coupled.** At maximum growth every strain can secrete zero α-KG — which is why the project proposes an addiction system.
4. **One deletion explains the design gap.** LitOpt out-produces the wet-lab build by +1.22 on glucose, and Δ*pta* (acetate overflow) accounts for +1.79 of it; Δ*ldhA* and Δ*poxB* contribute nothing in FBA. An earlier in-team analysis found *no* difference because it held the knockout background fixed and swapped only PPC/PYC — both results are correct and are reconciled in Figure 6.
5. **The plate endpoint is carbon-budget-limited.** In the 48-hour dFBA, wild-type (following true growth-optimal metabolism) secretes ~0, while both engineered strains clear the 8 mM target and converge near 22–23 mM. Carbon source sets the *timing* (glucose ~9 h, glycerol ~18 h), because glycerol supplies half the carbon per molecule.

## Install

Requires Python 3.10–3.13. The GLPK solver ships with COBRApy; nothing else to install.

```bash
git clone https://gitlab.igem.org/2026/software/houston-tx/flux-balance-analysis.git
cd flux-balance-analysis
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Run / reproduce the figures

One command regenerates all six wiki figures (PNG + SVG) into `figures/`:

```bash
python regenerate_model_figures.py
```

Or, equivalently, `make figures`. Runtime is roughly 3–6 minutes on a laptop (the dFBA and oxygen sweeps each solve several hundred LPs). The script prints every headline number before plotting, so you can verify the table above without opening an image:

```
O2 baseline: glucose=21.97, glycerol=13.53
Operating-point a-KG (50% O2, >=80% growth, mmol/gDW/h):
  glucose   WT=4.879  LitOpt=6.095  WetLab=4.879
  glycerol  WT=2.012  LitOpt=2.378  WetLab=2.012
PPC vs PYC (LitOpt): glucose 6.095/6.095  glycerol 2.378/2.378
```

The script is self-contained and reads `iHM1533.xml` from the repository root. Three switches at the top of the file control the modelling conventions:

| Setting | Default | Meaning |
|---|---|---|
| `O2_METHOD` | `"standard"` | O₂ baseline convention. `"pfba"` uses parsimonious FBA (lower baselines: 19.81 / 12.31). |
| `USE_PYC` | `False` | Model LitOpt's anaplerosis with native PPC. `True` swaps in heterologous PYC — identical result (see finding 2). |
| `GROWTH_FRACTION` | `0.80` | Growth floor for the operating point. |

## Repository layout

| Path | What it is |
|---|---|
| `regenerate_model_figures.py` | **Canonical figure generator.** Self-contained; produces all six wiki figures. |
| `figures/` | The six figures as published (PNG + SVG). |
| `iHM1533.xml`, `iHM1533.json` | The EcN genome-scale model (van 't Hof et al., 2022), SBML and JSON. |
| `docs/02_Model_wiki_page.md` | The Model wiki page text, with captions naming the files in `figures/`. |
| `docs/ANALYSIS.md` | The full 17-stage analysis narrative with embedded code and reasoning — the provenance for every design decision. |
| `UH iGEM 2026 Tutorial.ipynb` | Exploratory notebook walking through the same 17 stages interactively. |
| `scripts/verify_iHM1533_numbers.py` | Standalone check that the model loads and the headline fluxes reproduce. |
| `requirements.txt` | Pinned dependencies. |

### Which figures are canonical?

`regenerate_model_figures.py` and `figures/` are canonical and match the wiki. The notebook and `docs/ANALYSIS.md` are the **exploratory record**: they use the earlier *tier-based* framing (WT / Tier 1 / Tier 2 / Tier 2.5 / LitOpt) rather than the published *three-strain* framing (WT / LitOpt / wet-lab build), so the figures differ in presentation.

They do **not** disagree on the numbers. The analysis independently reports LitOpt = 6.0951 (Stage 14), PPC ≡ PYC at +0.0 % (Stage 13), CS activity = 8.0 % at 5 mM α-KG with 150 µM NADH and 20.0 % without (Stage 10), and an O₂ baseline of 21.97 (Stage 0) — all reproduced exactly by the script. Two cautions are documented in `docs/ANALYSIS.md`: the retired Stage 16 export `fig_pub_05_cs_feedback.png` is annotated "38 %", which does not reproduce from the constants the analysis itself states and is superseded by Stage 10; and **6.67** is the Tier 2.5 value, not LitOpt's, so any text quoting 6.67 for LitOpt is in error.

## Modelling caveats

FBA predicts stoichiometric ceilings, not kinetics. The model's maximum growth on glucose (~0.97 h⁻¹) is an upper bound a real culture will not reach, so absolute values run optimistic and are best read as relative comparisons. The O₂ baseline is degenerate (O₂ uptake at maximum growth is not unique), so "% of aerobic maximum" is only defined once a convention is fixed — see `O2_METHOD`. The citrate-synthase feedback model (Figure 4) is an illustrative kinetic layer on top of FBA, not a fitted result. Five such choices are documented in full in the script's header.

## Citation

If you use this work, please cite the repository (see `CITATION.cff`) and the underlying model:

> van 't Hof, M., Mohite, O. S., Monk, J. M., Weber, T., Palsson, B. O., & Sommer, M. O. A. (2022). High-quality genome-scale metabolic network reconstruction of probiotic bacterium *Escherichia coli* Nissle 1917. *BMC Bioinformatics*, 23, 566. DOI: [10.1186/s12859-022-05108-9](https://doi.org/10.1186/s12859-022-05108-9)

Key references for the design decisions:

- Chin, R. M., et al. (2014). The metabolite α-ketoglutarate extends lifespan by inhibiting ATP synthase and TOR. *Nature*, 510, 397–401. DOI: [10.1038/nature13264](https://doi.org/10.1038/nature13264)
- Chao, Y.-P., & Liao, J. C. (1993). Alteration of growth yield by overexpression of phosphoenolpyruvate carboxylase and phosphoenolpyruvate carboxykinase in *Escherichia coli*. *Appl. Environ. Microbiol.*, 59, 4261–4265. DOI: [10.1128/aem.59.12.4261-4265.1993](https://doi.org/10.1128/aem.59.12.4261-4265.1993)
- Gokarn, R. R., Eiteman, M. A., & Altman, E. (2000). Metabolic analysis of *Escherichia coli* in the presence and absence of the carboxylating enzymes. *Appl. Environ. Microbiol.*, 66, 1844–1850. DOI: [10.1128/aem.66.5.1844-1850.2000](https://doi.org/10.1128/aem.66.5.1844-1850.2000)
- Vemuri, G. N., Eiteman, M. A., & Altman, E. (2005). Physiological response of central metabolism in *E. coli* to deletion of pyruvate oxidase and introduction of heterologous pyruvate carboxylase. *Biotechnol. Bioeng.*, 90, 64–76. DOI: [10.1002/bit.20418](https://doi.org/10.1002/bit.20418)
- Stokell, D. J., et al. (2003). Probing the roles of key residues in the unique regulatory NADH binding site of type II citrate synthase of *E. coli*. *J. Biol. Chem.*, 278, 35435–35443. DOI: [10.1074/jbc.M302786200](https://doi.org/10.1074/jbc.M302786200)
- Eikmanns, B. J., et al. (1994). Nucleotide sequence, expression and transcriptional analysis of the *Corynebacterium glutamicum gltA* gene encoding citrate synthase. *Microbiology*, 140, 1817–1828. DOI: [10.1099/13500872-140-8-1817](https://doi.org/10.1099/13500872-140-8-1817)

A complete reference list is in `docs/ANALYSIS.md`.

## License

MIT — see [LICENSE](LICENSE). The model `iHM1533` is redistributed under the terms of its original publication; please cite van 't Hof et al. (2022).

## Contact

Austin Routt — UH iGEM 2026 Dry Lab lead.
