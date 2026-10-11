"""
Reproduce every alpha-KG number used on the AKGut Model page from iHM1533.
Requires: cobra (`pip install cobra`) and iHM1533.xml (in the project repo).
Run: python verify_iHM1533_numbers.py
"""
import cobra, warnings; warnings.filterwarnings("ignore")
from cobra import Reaction

from pathlib import Path
_MODEL = Path(__file__).resolve().parent.parent / "iHM1533.xml"  # repo root
m = cobra.io.read_sbml_model(str(_MODEL))
BIO="BIOMASS_EcN_iHM1533_core_59p80M"
GLC,GLYC,O2,AKG="EX_glc__D_e","EX_glyc_e","EX_o2_e","EX_akg_e"
LITOPT=("AKGDH","AKGDH2","LDH_D","PTAr","GLUDy","SUCDi","POX")   # dsucA x2, dldhA, dpta, dgdhA, dsdh, dpoxB
WET   =("GLUDy","AKGDH","AKGDH2")                               # dgdhA + sucA suppressed (no ldhA/pta)
TIER25=("AKGDH","AKGDH2","LDH_D","PTAr","GLUDy")

def setup(model, carbon, o2lb):
    model.reactions.get_by_id(GLC).lower_bound  = -10 if carbon=="glucose"  else 0
    model.reactions.get_by_id(GLYC).lower_bound = -10 if carbon=="glycerol" else 0
    model.reactions.get_by_id(AKG).lower_bound  = 0
    model.reactions.get_by_id(AKG).upper_bound  = 1000
    model.reactions.get_by_id(O2).lower_bound   = o2lb

def o2_baseline(carbon):
    with m:
        setup(m, carbon, -1000); m.objective=BIO
        return abs(m.optimize().fluxes[O2])

def max_akg_at_80(carbon, o2lb, kos=(), transport_correction=False):
    with m:
        setup(m, carbon, o2lb)
        for r in kos: m.reactions.get_by_id(r).bounds=(0,0)
        if transport_correction:                       # LitOpt: KgtP import-only + free exporter + glyoxylate off
            m.reactions.get_by_id("ICL").bounds=(0,0); m.reactions.get_by_id("MALS").bounds=(0,0)
            m.reactions.get_by_id("AKGt2rpp").lower_bound=0
            exp=Reaction("AKG_free_exporter"); exp.lower_bound=0; exp.upper_bound=1000
            exp.add_metabolites({m.metabolites.akg_c:-1, m.metabolites.akg_p:1}); m.add_reactions([exp])
        m.objective=BIO; mu=m.optimize().objective_value
        if mu is None or mu<0.01: return None
        m.reactions.get_by_id(BIO).lower_bound=0.8*mu
        m.objective=AKG; return m.optimize().objective_value

print(f"model: {len(m.reactions)} rxns / {len(m.metabolites)} mets / {len(m.genes)} genes  (expect 3143/2261/1533)")
gb = o2_baseline("glucose"); yb = o2_baseline("glycerol")
print(f"O2 baseline glucose={gb:.2f} (21.97) | glycerol growth check")
print("\n50% O2, >=80% growth, max EX_akg_e (mmol/gDW/h):")
print(f"  WT glucose                 {max_akg_at_80('glucose', -0.5*gb):.3f}   (README 4.88)")
print(f"  Tier2.5 glucose            {max_akg_at_80('glucose', -0.5*gb, TIER25):.3f}   (README 6.67)")
print(f"  LitOpt glucose (corrected) {max_akg_at_80('glucose', -0.5*gb, LITOPT, True):.3f}   (README 6.0951)")
print(f"  WT glycerol                {max_akg_at_80('glycerol',-0.5*yb):.3f}   (doc 2.01)")
print(f"  Tier2.5 glycerol           {max_akg_at_80('glycerol',-0.5*yb, TIER25):.3f}   (doc 2.61)")
print(f"  LitOpt glycerol (corrected){max_akg_at_80('glycerol',-0.5*yb, LITOPT, True):.3f}")
print(f"  WET-LAB build glycerol     {max_akg_at_80('glycerol',-0.5*yb, WET):.3f}   (== WT glycerol)")
