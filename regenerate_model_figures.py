#!/usr/bin/env python3
"""
regenerate_model_figures.py
===========================
Regenerates every alpha-ketoglutarate (a-KG) figure on the UH iGEM 2026 "Model"
wiki page. All numbers are computed from the genome-scale metabolic model
iHM1533 (van 't Hof et al. 2022, BMC Bioinformatics 23:566); place iHM1533.xml
next to this script and run it.

Three strains are compared on BOTH glucose and glycerol:
    WT      - wild-type E. coli Nissle 1917 (no modifications)
    LitOpt  - the computational design (glucose- and PYC-native)
    WetLab  - the wet-lab build "Proposed Strategy" (glycerol-native)

Requirements : cobra>=0.29, numpy, matplotlib   (LP solver: GLPK, bundled w/ cobra)
Usage        : python regenerate_model_figures.py
Outputs      : PNG + SVG for
    fig1_production_envelopes    growth vs a-KG trade-off      (glucose | glycerol)
    fig2_oxygen_sensitivity      a-KG vs O2 fraction           (glucose | glycerol)
    fig3_akg_operating_point     operating-point a-KG bars + PPC-vs-PYC design panel
    fig4_cs_feedback             citrate-synthase feedback kinetics (a-KG x NADH)
    fig5_dfba_plate              48 h plate: strains x carbons, vs the 8 mM target

================================================================================
STRAIN DEFINITIONS (mirror the wiki "Two-design comparison" table)
================================================================================
WT     : no modifications.

LitOpt : DsucA x2 (AKGDH + AKGDH2), DldhA (LDH_D), Dpta (PTAr), DgdhA (GLUDy),
         DsdhCDAB (SUCDi), DpoxB (POX); glyoxylate shunt held off (ICL = MALS = 0);
         native KgtP (AKGt2rpp) made import-only PLUS a dedicated free a-KG
         exporter; heterologous pyruvate carboxylase (PYC) for anaplerosis;
         feedback-insensitive citrate synthase (CggltA). PYC and CggltA are the
         design's two heterologous enzymes -- see notes (2) and (3).

WetLab : DgdhA (GLUDy) + ~80% sucA CRISPRi knockdown ONLY. It has NO dedicated
         a-KG exporter and does NOT delete ldhA / pta / poxB / sdh. gltA(R163L)
         and ppc overexpression and the optional DarcA are design elements whose
         effect is kinetic/regulatory -- see notes (2) and (4).

================================================================================
FIVE MODELLING CHOICES THAT MATTER (documented so the script is honest)
================================================================================
(1) THE O2 BASELINE IS DEGENERATE. O2 uptake at maximum growth is NOT unique in
    this model (by FVA, glucose spans 19.8-270 mmol/gDW/h). "% of aerobic maximum"
    is therefore only defined once a convention is fixed. We use the STANDARD-FBA
    solution (GLPK): 21.97 (glucose) / 13.53 (glycerol). Parsimonious FBA gives
    lower baselines (19.81 / 12.31); Shalaka's separate glycerol notebook used the
    pFBA glycerol value, which is why it reports ~1.87 rather than ~2.01 for the
    wet-lab build. Set O2_METHOD = "pfba" to switch conventions.

(2) CggltA AND gltA(R163L) ARE KINETIC, NOT STOICHIOMETRIC. Flux-balance analysis
    does not represent allosteric (NADH) or competitive (a-KG) inhibition, so
    relieving citrate-synthase feedback changes no FBA flux. Its effect is shown
    separately in Figure 4. Native citrate synthase is left unconstrained in every
    FBA calculation here -- that IS the "feedback-insensitive CS" assumption.

(3) PYC AND PPC ARE STOICHIOMETRICALLY EQUIVALENT IN FBA; PYC'S BENEFIT IS KINETIC
    (Figure 3B). LitOpt was designed around heterologous pyruvate carboxylase (PYC),
    which draws OAA from pyruvate + ATP instead of from PEP (native PPC). With a
    correctly mass- AND charge-balanced PYC reaction (it releases one H+), PYC and
    PPC give IDENTICAL a-KG in this model: 6.10 on glucose, 2.38 on glycerol. (An
    earlier version of this script used an unbalanced PYC missing that H+, which
    acted as a phantom proton source and spuriously inflated the value to 6.73 --
    corrected here.) The real reason to prefer PYC is KINETIC and FBA cannot see it:
    glucose uptake through the PTS consumes PEP, so overexpressing PPC competes with
    sugar uptake and lowers glucose consumption (Chao & Liao 1993), whereas the
    R. etli PYC replenishes OAA from pyruvate without touching PEP (Gokarn et al.
    2000), and DpoxB + PYC raises biomass yield ~42% and cuts acetate ~80%
    (Vemuri et al. 2005). LitOpt is therefore modelled with native PPC (USE_PYC =
    False) because the FBA number is identical either way; Figure 3B shows the
    equivalence and states the kinetic rationale.

(4) arcA IS A REGULATOR, NOT A METABOLIC REACTION, so the wet-lab build's optional
    DarcA cannot be represented in stoichiometric FBA. The wet-lab build is
    modelled without it. Its ~80% sucA CRISPRi knockdown is modelled as full
    suppression of AKGDH/AKGDH2: under a-KG maximisation a "<=20% residual capacity"
    cap collapses to zero flux, so an 80% knockdown and a full knockout give the
    identical a-KG ceiling (reproduces Shalaka's "+0.0000%"; see verify_iHM1533_numbers.py).

(5) THE dFBA PLATE ENDPOINT CANNOT DISCRIMINATE STRAINS -- BY DESIGN OF THE SYSTEM,
    NOT A BUG. Two facts, both reproduced here: (a) at maximum growth EVERY strain
    secretes ~0 a-KG (none is growth-coupled -- this is the rationale for the
    proposed addiction system); (b) once biomass saturates, the final plate [a-KG]
    is set by the plate CARBON BUDGET, so any operated producer converges to
    ~22-24 mM. In Figure 5 we therefore model WT by its true growth-optimal
    metabolism -- it re-oxidises a-KG through AKGDH and secretes ~0, the correct
    negative control -- and the engineered strains as operated producers. The
    result: WT makes nothing, both engineered strains clear the 8 mM Chin (2014)
    target, and CARBON SOURCE sets the timing (glucose ~9-10 h, glycerol ~18-19 h,
    because glycerol supplies half the carbon per molecule). The strain-vs-strain
    difference (LitOpt vs wet-lab) is the operating-point RATE and lives in
    Figures 1-3, NOT in the plate endpoint.
"""

import warnings
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cobra
from cobra import Reaction

# ------------------------------------------------------------------ config ----
# Paths are resolved relative to THIS file, so the script can be run from any
# working directory (`python regenerate_model_figures.py`, `make figures`, ...).
from pathlib import Path
REPO_ROOT  = Path(__file__).resolve().parent
MODEL_FILE = str(REPO_ROOT / "iHM1533.xml")   # the EcN model, at the repo root
OUT_DIR    = REPO_ROOT / "figures"            # all figures are written here
O2_METHOD  = "standard"    # "standard" (matches wiki) or "pfba"
USE_PYC    = False         # PPC and (balanced) PYC are stoichiometrically IDENTICAL (6.10); see note 3
GROWTH_FRACTION = 0.80     # operating point: >= 80% of a strain's own max growth
DPI        = 200

BIO   = "BIOMASS_EcN_iHM1533_core_59p80M"
EX_GLC, EX_GLYC, EX_O2, EX_AKG = "EX_glc__D_e", "EX_glyc_e", "EX_o2_e", "EX_akg_e"
KGTP  = "AKGt2rpp"

C_WT, C_LITOPT, C_WETLAB = "#8A8A8A", "#2CA02C", "#F07E1A"
STRAIN_COLOR = {"WT": C_WT, "LitOpt": C_LITOPT, "WetLab": C_WETLAB}
plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none"})


# =====================================================================
#  MODEL / STRAIN MACHINERY
# =====================================================================
def load_model():
    m = cobra.io.read_sbml_model(MODEL_FILE)
    assert (len(m.reactions), len(m.metabolites), len(m.genes)) == (3143, 2261, 1533), \
        "Unexpected model size -- is this iHM1533.xml?"
    return m


def set_medium(model, carbon, o2_lb):
    model.reactions.get_by_id(EX_GLC).lower_bound  = -10.0 if carbon == "glucose"  else 0.0
    model.reactions.get_by_id(EX_GLYC).lower_bound = -10.0 if carbon == "glycerol" else 0.0
    model.reactions.get_by_id(EX_AKG).lower_bound  = 0.0
    model.reactions.get_by_id(EX_AKG).upper_bound  = 1000.0
    model.reactions.get_by_id(EX_O2).lower_bound   = o2_lb


def apply_strain(model, strain, use_pyc=None):
    """Apply a strain's modifications IN PLACE (use inside `with model:`).
    use_pyc overrides the global USE_PYC (used by Figure 3B to compare PPC vs PYC)."""
    if use_pyc is None:
        use_pyc = USE_PYC
    if strain == "WT":
        return
    if strain == "LitOpt":
        for r in ("AKGDH", "AKGDH2", "LDH_D", "PTAr", "GLUDy", "SUCDi", "POX"):
            model.reactions.get_by_id(r).bounds = (0.0, 0.0)
        model.reactions.get_by_id("ICL").bounds  = (0.0, 0.0)     # glyoxylate shunt off
        model.reactions.get_by_id("MALS").bounds = (0.0, 0.0)
        model.reactions.get_by_id(KGTP).lower_bound = 0.0         # KgtP import-only
        exp = Reaction("AKG_free_exporter"); exp.bounds = (0.0, 1000.0)
        exp.add_metabolites({model.metabolites.akg_c: -1, model.metabolites.akg_p: 1})
        model.add_reactions([exp])
        if use_pyc:                                               # PYC replaces PPC (note 3)
            model.reactions.get_by_id("PPC").bounds = (0.0, 0.0)
            pyc = Reaction("PYC_hetero"); pyc.bounds = (0.0, 1000.0)
            pyc.add_metabolites({model.metabolites.pyr_c: -1, model.metabolites.hco3_c: -1,
                                 model.metabolites.atp_c: -1, model.metabolites.oaa_c: 1,
                                 model.metabolites.adp_c: 1,  model.metabolites.pi_c: 1,
                                 model.metabolites.h_c: 1})   # H+ product: mass- AND charge-balanced
            model.add_reactions([pyc])
        return
    if strain == "WetLab":
        model.reactions.get_by_id("GLUDy").bounds  = (0.0, 0.0)   # DgdhA
        model.reactions.get_by_id("AKGDH").bounds  = (0.0, 0.0)   # ~80% sucA CRISPRi
        model.reactions.get_by_id("AKGDH2").bounds = (0.0, 0.0)   #   (== full suppression, note 4)
        return
    raise ValueError(strain)


def o2_baseline(model, carbon):
    with model:
        set_medium(model, carbon, -1000.0)
        model.objective = BIO
        if O2_METHOD == "pfba":
            from cobra.flux_analysis import pfba
            return abs(pfba(model).fluxes[EX_O2])
        return abs(model.optimize().fluxes[EX_O2])


def max_akg_at_growth(model, strain, carbon, o2_lb, growth_fraction=GROWTH_FRACTION, use_pyc=None):
    """Operating-point value: hold growth >= growth_fraction * mu_max, maximise a-KG."""
    with model:
        set_medium(model, carbon, o2_lb)
        apply_strain(model, strain, use_pyc=use_pyc)
        model.objective = BIO
        mu = model.optimize().objective_value
        if mu is None or mu < 1e-2:
            return None
        model.reactions.get_by_id(BIO).lower_bound = growth_fraction * mu
        model.objective = EX_AKG
        v = model.optimize().objective_value
        return v if v is not None else None


def production_envelope(model, strain, carbon, o2_lb, n=26):
    with model:
        set_medium(model, carbon, o2_lb)
        apply_strain(model, strain)
        model.objective = BIO
        mu_max = model.optimize().objective_value
        gs, as_ = [], []
        for g in np.linspace(0.0, mu_max, n):
            model.reactions.get_by_id(BIO).bounds = (g, g)
            model.objective = EX_AKG
            v = model.optimize().objective_value
            gs.append(g); as_.append(v if v is not None else 0.0)
        return np.array(gs), np.array(as_)


# =====================================================================
#  dFBA (Static Optimization Approach; reproduces README Stage 15)
# =====================================================================
PLATE_VOL_L = 0.012
X_INIT_gDW  = 0.0001
V_MAX_UPT   = 10.0
KM_UPT_mM   = 0.5
DT_H        = 0.1
T_END_H     = 48.0
TARGET_mM   = 8.0


def run_dfba(model, strain, carbon, init_mM, x_max_gDW, o2_lb, produce=True):
    """SOA dFBA on a plate. produce=True -> operated producer (maximise a-KG at the
    operating point / in stationary phase); produce=False -> follow growth-optimal
    metabolism (used for WT, which re-oxidises a-KG and secretes ~0; see note 5).
    Returns (t, biomass_mg, carbon_mM, akg_mM)."""
    ex_c = EX_GLC if carbon == "glucose" else EX_GLYC
    n = int(T_END_H / DT_H) + 1
    t = np.linspace(0.0, T_END_H, n)
    X = np.zeros(n); C = np.zeros(n); A = np.zeros(n)
    X[0] = X_INIT_gDW
    C[0] = init_mM * PLATE_VOL_L
    with model:
        set_medium(model, carbon, o2_lb)
        apply_strain(model, strain)
        c_rxn = model.reactions.get_by_id(ex_c)
        b_rxn = model.reactions.get_by_id(BIO)
        c_lb0 = c_rxn.lower_bound
        for i in range(1, n):
            x, c, a = X[i-1], C[i-1], A[i-1]
            c_mM = c / PLATE_VOL_L
            if c_mM < 1e-4 or x < 1e-10:
                X[i], C[i], A[i] = x, max(c, 0.0), a
                continue
            v_upt = min(V_MAX_UPT * c_mM / (KM_UPT_mM + c_mM), c / (x * DT_H))
            growth_frac = max(1.0 - x / x_max_gDW, 0.0)
            mu = q_c = q_a = 0.0
            c_rxn.lower_bound = -v_upt
            if growth_frac > 0.01:                               # exponential phase
                b_rxn.bounds = (0.0, 1000.0)
                model.objective = BIO
                s1 = model.optimize()
                if s1.status == "optimal" and s1.objective_value > 1e-9:
                    ceil = s1.objective_value * growth_frac
                    if produce:                                   # operated producer
                        b_rxn.bounds = (0.80 * ceil, ceil)
                        model.objective = EX_AKG
                        s2 = model.optimize()
                    else:                                         # WT: grow, natural a-KG (~0)
                        b_rxn.bounds = (ceil, ceil)
                        model.objective = BIO
                        s2 = model.optimize()
                    if s2.status == "optimal":
                        mu, q_c = s2.fluxes[BIO], s2.fluxes[ex_c]
                        q_a = s2.fluxes.get(EX_AKG, 0.0)
                    else:
                        q_c = -v_upt
                else:
                    q_c = -v_upt
            else:                                                # stationary phase
                b_rxn.bounds = (0.0, 0.0)
                model.objective = EX_AKG if produce else BIO
                ss = model.optimize()
                if ss.status == "optimal":
                    q_c = ss.fluxes[ex_c]
                    q_a = ss.fluxes.get(EX_AKG, 0.0) if produce else 0.0
                else:
                    q_c = -v_upt
            c_rxn.lower_bound = c_lb0
            b_rxn.bounds = (0.0, 1000.0)
            X[i] = min(x + mu * x * DT_H, x_max_gDW)
            C[i] = max(c + q_c * x * DT_H, 0.0)
            A[i] = max(a + q_a * x * DT_H, 0.0)
    return t, X * 1000.0, C / PLATE_VOL_L, A / PLATE_VOL_L


# =====================================================================
#  CITRATE-SYNTHASE FEEDBACK KINETICS (reproduces README Stage 10)
# =====================================================================
KM_OAA, OAA_mM, KI_AKG, KI_NADH = 0.020, 0.005, 1.0, 0.100


def cs_rel_activity(akg_mM, nadh_mM):
    comp = (KM_OAA + OAA_mM) / (KM_OAA * (1.0 + akg_mM / KI_AKG) + OAA_mM)
    allo = 1.0 / (1.0 + nadh_mM / KI_NADH)
    return 100.0 * comp * allo


# =====================================================================
#  FIGURES
# =====================================================================
def _save(fig, name):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stem = OUT_DIR / name
    fig.savefig(str(stem) + ".png", dpi=DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(str(stem) + ".svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote figures/{name}.png / .svg")


def fig1_production_envelopes(model, base):
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6))
    for ax, carbon in zip(axes, ("glucose", "glycerol")):
        o2 = -0.5 * base[carbon]
        for strain in ("WT", "LitOpt", "WetLab"):
            g, a = production_envelope(model, strain, carbon, o2)
            ax.plot(g, a, color=STRAIN_COLOR[strain], lw=2.4, label=strain, zorder=3)
            pt = max_akg_at_growth(model, strain, carbon, o2)
            with model:
                set_medium(model, carbon, o2); apply_strain(model, strain)
                model.objective = BIO; mu = model.optimize().objective_value
            if pt is not None:
                ax.plot(0.8 * mu, pt, "*", color=STRAIN_COLOR[strain], ms=15,
                        markeredgecolor="white", markeredgewidth=0.8, zorder=5)
        ax.set_title(f"{carbon.capitalize()}  (50% O\u2082)", fontsize=13, fontweight="bold")
        ax.set_xlabel("Growth rate  (h\u207b\u00b9)", fontsize=12)
        ax.set_ylabel("Max \u03b1-KG secretion  (mmol\u00b7gDW\u207b\u00b9\u00b7h\u207b\u00b9)", fontsize=12)
        ax.legend(fontsize=10.5, loc="upper right")
        for s in ("top", "right"): ax.spines[s].set_visible(False)
    fig.suptitle("Production envelopes: growth vs \u03b1-KG (WT / LitOpt / wet-lab build)",
                 fontsize=15, fontweight="bold", y=1.02)
    fig.text(0.5, -0.03, "\u2605 = 80% growth operating point.",
             ha="center", fontsize=9.5, style="italic", color="#555555")
    fig.tight_layout()
    _save(fig, "fig1_production_envelopes")


def fig2_oxygen_sensitivity(model, base):
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6))
    fracs = np.linspace(0.0, 1.0, 21)
    for ax, carbon in zip(axes, ("glucose", "glycerol")):
        for strain in ("WT", "LitOpt", "WetLab"):
            ys = [max_akg_at_growth(model, strain, carbon, -f * base[carbon]) for f in fracs]
            ys = [v if v is not None else np.nan for v in ys]
            ax.plot(fracs * 100, ys, color=STRAIN_COLOR[strain], lw=2.4, marker="o", ms=3.5,
                    label=strain, zorder=3)
        ax.axvline(50, color="#cc4444", ls="--", lw=1.3, zorder=1)
        ax.text(51, ax.get_ylim()[1] * 0.30, "50% (operating point)", color="#cc4444",
                fontsize=8.5, rotation=90, va="bottom", ha="left")
        ax.set_title(f"{carbon.capitalize()}", fontsize=13, fontweight="bold")
        ax.set_xlabel("O\u2082 uptake  (% of aerobic maximum)", fontsize=12)
        ax.set_ylabel("Max \u03b1-KG at \u226580% growth  (mmol\u00b7gDW\u207b\u00b9\u00b7h\u207b\u00b9)", fontsize=12)
        ax.legend(fontsize=10.5, loc="upper left")
        for s in ("top", "right"): ax.spines[s].set_visible(False)
    fig.suptitle("Oxygen sensitivity of \u03b1-KG production", fontsize=15, fontweight="bold", y=1.02)
    fig.text(0.5, -0.03, "Non-monotonic 'hump': production peaks near 50% O\u2082 (microaerobic).",
             ha="center", fontsize=9.5,
             style="italic", color="#555555")
    fig.tight_layout()
    _save(fig, "fig2_oxygen_sensitivity")


def fig3_akg_operating_point(model, base):
    strains = ("WT", "LitOpt", "WetLab")
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(14.0, 6.2))

    # Panel A: WT / LitOpt / WetLab at the operating point, both carbons
    data = {c: [max_akg_at_growth(model, s, c, -0.5 * base[c]) for s in strains]
            for c in ("glucose", "glycerol")}
    width = 0.24; xg = np.arange(2) * 0.9
    for k, s in enumerate(strains):
        vals = [data["glucose"][k], data["glycerol"][k]]
        bars = axA.bar(xg + (k-1)*width, vals, width=width*0.95, color=STRAIN_COLOR[s],
                       edgecolor="white", label=s, zorder=3)
        for b, v in zip(bars, vals):
            axA.text(b.get_x()+b.get_width()/2, v+0.08, f"{v:.2f}", ha="center", va="bottom",
                     fontsize=10.5, fontweight="bold", color="#333")
    axA.set_xticks(xg); axA.set_xticklabels(["Glucose", "Glycerol"], fontsize=12.5)
    axA.set_ylabel("Max \u03b1-KG at 50% O\u2082, \u226580% growth\n(mmol\u00b7gDW\u207b\u00b9\u00b7h\u207b\u00b9)", fontsize=12)
    axA.set_title("A. Strain comparison", fontsize=13, fontweight="bold")
    axA.legend(fontsize=10.5, loc="upper right")
    axA.set_ylim(0, max(max(data["glucose"]), max(data["glycerol"])) * 1.22)
    for s in ("top", "right"): axA.spines[s].set_visible(False)

    # Panel B: PPC vs PYC for LitOpt. With a balanced PYC these are IDENTICAL in FBA;
    # PYC's real advantage is kinetic (PTS protection) and cannot appear in FBA (note 3).
    ppc = [max_akg_at_growth(model, "LitOpt", c, -0.5 * base[c], use_pyc=False) for c in ("glucose", "glycerol")]
    pyc = [max_akg_at_growth(model, "LitOpt", c, -0.5 * base[c], use_pyc=True)  for c in ("glucose", "glycerol")]
    xg2 = np.arange(2) * 0.85; w2 = 0.32
    b1 = axB.bar(xg2 - w2/2, ppc, width=w2, color="#9ecae1", edgecolor="white", label="native PPC", zorder=3)
    b2 = axB.bar(xg2 + w2/2, pyc, width=w2, color="#2CA02C", edgecolor="white", label="heterologous PYC", zorder=3)
    for bars in (b1, b2):
        for b in bars:
            axB.text(b.get_x()+b.get_width()/2, b.get_height()+0.06, f"{b.get_height():.2f}",
                     ha="center", va="bottom", fontsize=10, fontweight="bold", color="#333")
    for x in xg2:                              # "=" between the identical bars
        axB.text(x, max(ppc[0], pyc[0]) * 0.5, "=", ha="center", va="center",
                 fontsize=20, fontweight="bold", color="#1a7a1a")
    axB.text(0.5, 0.94, "Identical in FBA \u2014 PYC's benefit is KINETIC (PTS protection),\n"
             "which FBA cannot represent: PPC overexpression competes with glucose\n"
             "uptake for PEP (Chao & Liao 1993); R. etli PYC spares PEP (Gokarn 2000);\n"
             "\u0394poxB + PYC \u2192 +42% biomass, \u221280% acetate (Vemuri 2005).",
             transform=axB.transAxes, ha="center", va="top", fontsize=8.2, color="#333",
             bbox=dict(boxstyle="round", fc="#f2f8f2", ec="#bcd8bc"))
    axB.set_xticks(xg2); axB.set_xticklabels(["Glucose", "Glycerol"], fontsize=12.5)
    axB.set_ylabel("LitOpt max \u03b1-KG  (mmol\u00b7gDW\u207b\u00b9\u00b7h\u207b\u00b9)", fontsize=12)
    axB.set_title("B. Anaplerosis: PPC vs PYC (stoichiometrically equivalent)", fontsize=12, fontweight="bold")
    axB.legend(fontsize=10, loc="center right", bbox_to_anchor=(1.0, 0.53))
    axB.set_ylim(0, max(pyc) * 1.7)
    for s in ("top", "right"): axB.spines[s].set_visible(False)

    fig.suptitle("\u03b1-KG at the operating point: carbon source, design, and anaplerosis",
                 fontsize=15, fontweight="bold", y=1.03)
    fig.text(0.5, -0.02, "The wet-lab build \u2248 WT on both carbons (it omits the overflow knockouts + exporter). "
             "In FBA, PPC and PYC give the same yield; PYC is chosen for its kinetic PTS benefit.",
             ha="center", fontsize=9.2, style="italic", color="#555555")
    fig.tight_layout()
    _save(fig, "fig3_akg_operating_point")


def fig4_cs_feedback():
    akg = np.linspace(0, 10, 400)
    series = [(0.0, "#9AA0A6", "NADH = 0 (\u03b1-KG only)", 1.8),
              (0.080, "#3E8FB8", "NADH = 80 \u00b5M", 1.8),
              (0.150, "#D6336C", "NADH = 150 \u00b5M  (producing strain)", 3.2),
              (0.300, "#5B3A8C", "NADH = 300 \u00b5M", 1.8)]
    fig, ax = plt.subplots(figsize=(10.5, 7.0))
    for nadh, c, lbl, lw in series:
        ax.plot(akg, [cs_rel_activity(x, nadh) for x in akg], color=c, lw=lw, label=lbl, zorder=3)
    ax.axvspan(1, 10, color="#EAF4EA", alpha=0.6, zorder=0)
    ax.text(9.7, 96, "expected producer\n[\u03b1-KG]\u1d62 range", ha="right", va="top", fontsize=9, color="#4a7a4a")
    ax.axvline(0.5, color="#3E8FB8", ls="--", lw=1.3); ax.text(0.6, 5, "WT\n\u22480.5 mM", fontsize=9, color="#3E8FB8", va="bottom")
    ax.axvline(5.0, color="#D6336C", ls="--", lw=1.3)
    y8 = cs_rel_activity(5.0, 0.150)
    ax.plot(5.0, y8, "o", color="#D6336C", ms=8, zorder=5)
    ax.annotate(f"\u2248{y8:.0f}% activity\nat 5 mM \u03b1-KG,\n150 \u00b5M NADH", xy=(5.0, y8),
                xytext=(6.4, 26), fontsize=11.5, fontweight="bold", color="#D6336C",
                arrowprops=dict(arrowstyle="->", color="#D6336C", lw=1.8))
    ax.annotate("NADH alone caps activity\nat 40% (150 \u00b5M) before\nany \u03b1-KG accumulates",
                xy=(0.05, 40), xytext=(2.0, 64), fontsize=10, color="#8a1e46",
                arrowprops=dict(arrowstyle="->", color="#8a1e46", lw=1.4))
    ax.set_xlim(0, 10); ax.set_ylim(0, 102)
    ax.set_xlabel("Intracellular [\u03b1-KG]  (mM)", fontsize=13)
    ax.set_ylabel("Relative citrate-synthase activity  (%)", fontsize=13)
    ax.set_title("Citrate synthase feedback: \u03b1-KG competition \u00d7 NADH allosteric inhibition",
                 fontsize=14.5, fontweight="bold", pad=26)
    ax.text(0.5, 1.03, "Illustrative kinetic model \u2014 E. coli type II citrate synthase is allosterically "
            "inhibited by NADH (Stokell et al., 2003); motivates the feedback-insensitive CggltA swap",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=9.3, style="italic", color="#555555")
    ax.legend(loc="upper right", fontsize=10.5, framealpha=0.95)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    _save(fig, "fig4_cs_feedback")


def fig5_dfba_plate(model, base):
    """Strain x carbon plate simulation (see note 5). WT follows growth-optimal
    metabolism (a-KG ~ 0, the negative control); engineered strains are operated
    producers. Rows = carbon source; columns = biomass, substrate, a-KG."""
    carbons = [("glucose", 22.2), ("glycerol", 43.4)]     # 0.4% (w/v)
    fig, axes = plt.subplots(2, 3, figsize=(16.0, 9.0))
    col_titles = ["Lawn growth", "Substrate depletion", "\u03b1-KG accumulation"]
    for r, (carbon, init_mM) in enumerate(carbons):
        o2 = -0.5 * base[carbon]
        for strain in ("WT", "LitOpt", "WetLab"):
            produce = (strain != "WT")
            t, X, C, A = run_dfba(model, strain, carbon, init_mM, 0.006, o2, produce=produce)
            lbl = strain + (" (control, \u22480)" if strain == "WT" else "")
            axes[r, 0].plot(t, X, color=STRAIN_COLOR[strain], lw=2.4, label=lbl)
            axes[r, 1].plot(t, C, color=STRAIN_COLOR[strain], lw=2.4, label=strain)
            axes[r, 2].plot(t, A, color=STRAIN_COLOR[strain], lw=2.4,
                            label=f"{strain}: {A[-1]:.0f} mM" if produce else "WT: \u22480 mM")
        axes[r, 2].axhline(TARGET_mM, color="#cc2222", ls="--", lw=1.8)
        axes[r, 2].text(24, TARGET_mM + 0.7, "8 mM target", color="#cc2222", fontsize=9)
        axes[r, 0].set_ylabel(f"{carbon.capitalize()}\n\nBiomass (mg)", fontsize=12, fontweight="bold")
        axes[r, 1].set_ylabel("[substrate] (mM)", fontsize=11)
        axes[r, 2].set_ylabel("[\u03b1-KG] (mM)", fontsize=11)
        for k in range(3):
            axes[r, k].set_xlabel("Time (h)", fontsize=11)
            axes[r, k].legend(fontsize=8.5, loc="best")
            for s in ("top", "right"): axes[r, k].spines[s].set_visible(False)
            if r == 0:
                axes[r, k].set_title(col_titles[k], fontsize=13, fontweight="bold")
    fig.suptitle("dFBA plate simulation (SOA): strains \u00d7 carbon sources on 12 mL modified NGM + 0.4% carbon",
                 fontsize=15, fontweight="bold", y=1.005)
    fig.text(0.5, -0.015,
             "WT (control) re-oxidises \u03b1-KG and secretes ~0; both engineered strains clear 8 mM. The plate ENDPOINT is "
             "carbon-budget-limited, so LitOpt and the wet-lab build converge (~22-23 mM) \u2014 their design difference is the "
             "operating-point RATE (Figures 1-3). Carbon source sets the TIMING (glucose ~9-10 h, glycerol ~18-19 h: half the "
             "carbon per molecule).",
             ha="center", fontsize=9.0, style="italic", color="#555555")
    fig.tight_layout()
    _save(fig, "fig5_dfba_plate")


# =====================================================================
#  MAIN
# =====================================================================
def main():
    print("Loading", MODEL_FILE, "...")
    model = load_model()
    print(f"  model OK: {len(model.reactions)} rxns / {len(model.metabolites)} mets / {len(model.genes)} genes")
    print(f"  config: O2_METHOD={O2_METHOD}, USE_PYC={USE_PYC}")

    base = {c: o2_baseline(model, c) for c in ("glucose", "glycerol")}
    print(f"O2 baseline: glucose={base['glucose']:.2f}, glycerol={base['glycerol']:.2f}")

    print("\nOperating-point a-KG (50% O2, >=80% growth, mmol/gDW/h):")
    for c in ("glucose", "glycerol"):
        row = {s: max_akg_at_growth(model, s, c, -0.5 * base[c]) for s in ("WT", "LitOpt", "WetLab")}
        print(f"  {c:8s}  WT={row['WT']:.3f}  LitOpt={row['LitOpt']:.3f}  WetLab={row['WetLab']:.3f}")
    print("PPC vs PYC (LitOpt): glucose {:.3f}/{:.3f}  glycerol {:.3f}/{:.3f}".format(
        max_akg_at_growth(model, "LitOpt", "glucose", -0.5*base["glucose"], use_pyc=False),
        max_akg_at_growth(model, "LitOpt", "glucose", -0.5*base["glucose"], use_pyc=True),
        max_akg_at_growth(model, "LitOpt", "glycerol", -0.5*base["glycerol"], use_pyc=False),
        max_akg_at_growth(model, "LitOpt", "glycerol", -0.5*base["glycerol"], use_pyc=True)))

    print("\nGenerating figures ...")
    fig1_production_envelopes(model, base)
    fig2_oxygen_sensitivity(model, base)
    fig3_akg_operating_point(model, base)
    fig4_cs_feedback()
    fig5_dfba_plate(model, base)
    print("\nDone. 5 figures written (PNG + SVG).")


if __name__ == "__main__":
    main()
