# Model

**Track:** Dry Lab · **Owner:** Austin Routt · **Wiki:** https://2026.igem.wiki/houston-tx/model · **Repository:** https://github.com/Austin-Routt/iGEM-2026

## The Central Question

Every cell growing on sugar makes α-ketoglutarate (α-KG) as a routine step of the tricarboxylic acid (TCA) cycle — the loop cells use to extract energy from food — then immediately consumes it again. We want to interrupt that consumption so *Escherichia coli* Nissle 1917 (EcN), our probiotic chassis, leaks α-KG onto a worm plate while still growing into a lawn. Before touching a pipette, we asked a computer model one question: what genetic edits force α-KG out of the cell without stalling growth? Two follow-ups shaped everything, because the wet lab chose a different recipe than our first design — how does the answer change on glycerol instead of glucose, and when the main α-KG-consuming gene (*sucA*) is dialed down ~80% rather than deleted?

## The FBA Approach

We answered with flux balance analysis (FBA), which predicts how carbon flows through metabolism from stoichiometry alone, without hard-to-measure reaction speeds. It represents the cell as a genome-scale model — a full inventory of every reaction the organism can run — and finds the steady-state flux pattern that maximizes a chosen objective, usually growth, within limits such as nutrient uptake. Its blind spots matter: FBA sees stoichiometry but not enzyme kinetics or regulation, so its numbers are ceilings, not guarantees.

We used the published EcN model iHM1533 (van 't Hof et al., 2022) — 3,143 reactions, 2,261 metabolites, 1,533 genes — in COBRApy. A *production envelope* holds growth above a set fraction and reports the most α-KG the cell can still make; *flux variability analysis* (FVA) asks whether α-KG is ever forced across equally optimal solutions. One calibration note keeps the page honest: the model's top predicted growth on glucose is ~0.97 h⁻¹, a stoichiometric ceiling real EcN never reaches, so absolute values run optimistic and are read as relative.

## The Analysis

The seventeen-stage analysis is one investigation. It opens with a diagnosis: the cell drains α-KG two ways at once — burning it through α-KG dehydrogenase (the *sucA* enzyme, 4.87 flux units) and converting it to glutamate through glutamate dehydrogenase (*gdhA*, 7.83 units). The larger, unexpected leak was glutamate dehydrogenase, which moved the target off *sucA* alone.

From there we built a ladder of strains. Tier 1 knocks out both α-KG-dehydrogenase paralogs; Tier 2 also seals the lactate and acetate overflow leaks (*ldhA*, *pta*); Tier 2.5 removes the biggest drain, *gdhA*; Tier 3 blocks glutamate synthesis entirely and kills the cell.

Oxygen gave the first surprise. Sweeping it from zero to full reveals a hump — α-KG peaks not at full oxygen but near half, the microaerobic regime, as partial limitation backs carbon into the TCA cycle (Figure 2). At this 50% operating point, Tier 2.5 tops the ladder at 6.67 units while still growing at 0.556 h⁻¹.

Then the catch: FVA showed no tier is growth-coupled. At maximum growth the cell can make zero α-KG, so production always costs growth, and faster non-producers would outcompete producers on a plate. A purely metabolic design is not enough.

Two enzyme-level problems followed. Citrate synthase, the gateway into the TCA cycle, is throttled by its own products; under production an illustrative kinetic model puts it near 8% activity (Figure 4). The fix is a feedback-insensitive citrate synthase from *Corynebacterium glutamicum* (CggltA), which is naturally NADH-insensitive (Eikmanns et al., 1994). Separately, the native α-KG transporter would re-import secreted product, so we run it import-only and add a dedicated exporter.

The remaining stages hardened the design. To counter growth-coupling we designed an addiction system — an essential gene under α-KG-responsive control — which works in the microaerobic gut but not the aerobic plate, so plate cultures need operational safeguards instead. The glyoxylate shunt stays off on glucose unaided, and the two anaplerotic (TCA-refilling) enzymes, PEP carboxylase (PPC) and pyruvate carboxylase (PYC), give identical FBA output — their difference is kinetic, discussed below.

## Key Findings

Three results carry the page. Tier 2.5 at 50% oxygen is the production sweet spot (6.67 units) — enough edits to redirect carbon, enough oxygen limitation to hold it, without starving growth. The feedback-insensitive citrate-synthase swap is not optional: without it, product throttles the enzyme feeding α-KG synthesis to ~8%. And because no design is growth-coupled, an addiction system is required to hold production over time, useful in the gut rather than on the plate.

## Predictions for the Wet Lab

The model makes testable predictions. Knocking down α-KG dehydrogenase (Δ*sucA*, both paralogs) is necessary but not sufficient; blocking lactate and acetate overflow (Δ*ldhA*, Δ*pta*) should raise measurable α-KG; removing glutamate dehydrogenase (Δ*gdhA*) clears the largest drain while an alternate route (GOGAT) keeps the cell prototrophic. After Δ*pta*, deleting pyruvate oxidase (Δ*poxB*) should close a secondary acetate leak — Vemuri et al. (2005) report ~42% higher biomass and ~80% less acetate in a *poxB*-*pyc*⁺ strain.

Attenuating rather than deleting *sucA* is the safer route, since succinyl-CoA supply limits growth in TCA-deleted *E. coli* (Zhou et al., 2024), which supports the wet lab's CRISPRi plan. The lab should measure growth versus titer, supernatant acetate and lactate, and α-KG at 20%, 50%, and 100% oxygen to test the predicted hump.

## The LitOpt Baseline

LitOpt (Literature-Optimized) is the complete glucose design: Δ*sucA*×2, Δ*ldhA*, Δ*pta*, Δ*gdhA*, Δ*sdhCDAB*, and Δ*poxB*, plus the CggltA swap, the import-only transporter with a dedicated exporter, the glyoxylate shunt held off, and heterologous PYC. It holds ~6.1 units of α-KG at the 50% oxygen operating point (Figures 1 and 3), below the uncorrected 6.67 because the import-only transport correction is deliberately conservative. Modeling it with PYC instead of PPC gives the same 6.1 — the two are stoichiometrically identical in FBA.

Fed into a 48-hour plate simulation, LitOpt reaches ~22.6 mM α-KG, about 2.8× the 8 mM at which Chin et al. (2014) saw maximal *C. elegans* lifespan extension (Figure 5). Two caveats attach: the plate endpoint is set by the carbon budget, so it is a proof-of-concept ceiling any producer shares, not a strain-specific number; and the strain is not growth-coupled. LitOpt's glucose-plus-PYC logic mirrors a published high-titer strain (Chen et al., 2020).

## From Glucose to Glycerol

The strain at the bench differs from LitOpt: it grows on glycerol, dials *sucA* down ~80% by CRISPRi (tunable gene repression) rather than deleting it, overexpresses a feedback-relieved citrate synthase (R163L) and *ppc* from a plasmid, and omits several knockouts — mirroring a published crude-glycerol strain (Chiang et al., 2025).

On glycerol the modeled yields run lower across every strain (Tier 2.5 near 2.6 versus 6.67), but most of that gap is a modeling convention, not a property of glycerol. Both carbons are supplied at the same molar uptake, and glycerol (three carbons) delivers half the carbon of glucose (six). At equal carbon, wild-type glycerol reaches 82% of the glucose value; the small remainder is glycerol's higher oxygen demand.

The PPC-versus-PYC choice is kinetic, not stoichiometric — which is why FBA rates them equal. Glucose import through the phosphotransferase system (PTS) spends one PEP per sugar (Gokarn et al., 2000), so overexpressing PPC starves uptake — glucose consumption falls up to ~40% at full induction (Chao & Liao, 1993) — whereas PYC draws on pyruvate and experimentally out-consumes PPC. Glycerol bypasses the PTS, so PPC is fine there.

Two model findings support the remaining wet-lab choices. An 80% *sucA* knockdown and a full knockout give identical α-KG in the model, backing tunable CRISPRi as the lower-risk route.

The feedback-relieved citrate synthase is also warranted, with a caveat: R163L weakens the enzyme's NADH inhibition more than 100-fold (Stokell et al., 2003) but not its α-KG inhibition, so the model's assumption that R163L relieves α-KG feedback is untested — unlike the established insensitivity of the *C. glutamicum* enzyme. Full design reasoning lives on our **Engineering** page; containment continues on **Safety**.

## Why two of our analyses disagreed

An earlier in-team analysis compared the two designs and found **no difference**; this page shows one. Both are correct — they ask different questions — and the reconciliation is itself a result (Figure 6).

That analysis held the knockout background fixed (a Tier 2.5 background already containing Δ*ldhA*, Δ*pta* and Δ*sdhCDAB*) and swapped only the anaplerotic enzyme, PPC against PYC. Both arms give 6.0951, which we reproduce exactly. That conclusion is sound: PPC and PYC are stoichiometrically equivalent.

This page instead takes each genotype as specified. The wet-lab build contains no Δ*ldhA*, Δ*pta*, Δ*sdhCDAB*, Δ*poxB* and no dedicated exporter, so the earlier test had effectively granted it LitOpt's knockouts before the comparison began.

Decomposing the gap edit by edit localises almost all of it to **one deletion, Δ*pta***, which blocks acetate overflow: it supplies +1.79 of the +1.22 net difference on glucose, the transport correction giving back −0.65. Δ*ldhA* and Δ*poxB* contribute nothing in FBA, Δ*sdhCDAB* almost nothing. The prediction is concrete: **adding Δ*pta* should capture most of LitOpt's advantage** — a single deletion in a strain already doing CRISPR knockouts.

## Figures

All figures are original (team-produced), regenerated from iHM1533 by `regenerate_model_figures.py`; every plotted value is computed from the model. Filenames are the files in `figures/`. Each compares wild-type (WT), LitOpt and the wet-lab build on glucose and glycerol.

- **Figure 1 — Production envelopes** (`figures/fig1_production_envelopes.png`). Growth versus maximum α-KG at 50% oxygen, 80% operating point marked. LitOpt sits highest; the wet-lab build tracks WT.
- **Figure 2 — Oxygen sensitivity** (`figures/fig2_oxygen_sensitivity.png`). Maximum α-KG at ≥80% growth across oxygen levels; the non-monotonic peak near 50% defines the operating point.
- **Figure 3 — Operating point and anaplerosis** (`figures/fig3_akg_operating_point.png`). (A) operating-point α-KG for all three strains; (B) PPC versus PYC for LitOpt — identical in FBA, annotated with the kinetic rationale FBA cannot capture.
- **Figure 4 — Citrate-synthase feedback** (`figures/fig4_cs_feedback.png`). Illustrative kinetic model; activity falls to ~8% under producing conditions, motivating the CggltA swap.
- **Figure 5 — dFBA plate simulation** (`figures/fig5_dfba_plate.png`). Biomass, substrate and α-KG over 48 hours. WT (control) makes ~0; both engineered strains clear 8 mM; the carbon-limited endpoint converges while carbon source sets timing (glucose ~9 h, glycerol ~18 h).
- **Figure 6 — Where the two designs differ** (`figures/fig6_strain_gap_decomposition.png`). (A) Adding LitOpt's edits to the wet-lab build one at a time: Δ*pta* alone gives +1.79 (glucose) and +0.60 (glycerol); Δ*ldhA* and Δ*poxB* add nothing. (B) Swapping only the anaplerotic enzyme gives no difference; comparing the genotypes as specified gives +1.22.

## References

Chao, Y.-P., & Liao, J. C. (1993). Alteration of growth yield by overexpression of phosphoenolpyruvate carboxylase and phosphoenolpyruvate carboxykinase in *Escherichia coli*. *Applied and Environmental Microbiology*, 59, 4261–4265. DOI: 10.1128/aem.59.12.4261-4265.1993

Chen, X., Dong, X., Liu, J., Luo, Q., & Liu, L. (2020). Pathway engineering of *Escherichia coli* for α-ketoglutaric acid production. *Biotechnology and Bioengineering*, 117, 2791–2801. DOI: 10.1002/bit.27456

Chiang, C.-J., Lin, J.-H., Ta, D.-T., Luu, N. L., Sun, Y.-C., & Chao, Y.-P. (2025). Metabolic engineering of *Escherichia coli* for overproduction of alpha-ketoglutarate using crude glycerol. *Journal of Agricultural and Food Chemistry*, 73, 18346–18352. DOI: 10.1021/acs.jafc.5c04670

Chin, R. M., Fu, X., Pai, M. Y., et al. (2014). The metabolite α-ketoglutarate extends lifespan by inhibiting ATP synthase and TOR. *Nature*, 510, 397–401. DOI: 10.1038/nature13264

Eikmanns, B. J., Thum-Schmitz, N., Eggeling, L., Lüdtke, K.-U., & Sahm, H. (1994). Nucleotide sequence, expression and transcriptional analysis of the *Corynebacterium glutamicum gltA* gene encoding citrate synthase. *Microbiology*, 140, 1817–1828. DOI: 10.1099/13500872-140-8-1817

Gokarn, R. R., Eiteman, M. A., & Altman, E. (2000). Metabolic analysis of *Escherichia coli* in the presence and absence of the carboxylating enzymes phosphoenolpyruvate carboxylase and pyruvate carboxylase. *Applied and Environmental Microbiology*, 66, 1844–1850. DOI: 10.1128/aem.66.5.1844-1850.2000

Stokell, D. J., Donald, L. J., Maurus, R., et al. (2003). Probing the roles of key residues in the unique regulatory NADH binding site of type II citrate synthase of *Escherichia coli*. *Journal of Biological Chemistry*, 278, 35435–35443. DOI: 10.1074/jbc.M302786200

van 't Hof, M., Mohite, O. S., Monk, J. M., Weber, T., Palsson, B. O., & Sommer, M. O. A. (2022). High-quality genome-scale metabolic network reconstruction of probiotic bacterium *Escherichia coli* Nissle 1917. *BMC Bioinformatics*, 23, 566. DOI: 10.1186/s12859-022-05108-9

Vemuri, G. N., Eiteman, M. A., & Altman, E. (2005). Physiological response of central metabolism in *Escherichia coli* to deletion of pyruvate oxidase and introduction of heterologous pyruvate carboxylase. *Biotechnology and Bioengineering*, 90, 64–76. DOI: 10.1002/bit.20418

Zhou, H., et al. (2024). A citric acid cycle-deficient *Escherichia coli* as an efficient chassis for aerobic fermentations. *Nature Communications*, 15, 2372. DOI: 10.1038/s41467-024-46655-4
