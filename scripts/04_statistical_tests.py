"""04 - Statistical comparison of formulations.

Inputs : results/tables/adsorption_by_replica.csv, adsorption_per_molecule.csv
Outputs: results/tables/statistical_tests.csv

Two levels of comparison are reported: one-way ANOVA on replica means, which
treats the independent simulation as the unit of replication (n = 3), and
non-parametric tests on individual molecules (n = 24 per condition), which have
more power but treat molecules within a replica as independent.
"""
import sys
from pathlib import Path

import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from algchi_screen.config import SYSTEMS, SYSTEM_LABEL, TABLES

PAIRS = [
    ("1-2", "ATC", "1-1", "ATC"),
    ("1-2", "MTC", "1-1", "MTC"),
    ("1-2", "CPZ", "1-1", "CPZ"),
    ("1-2", "ATC", "2-1", "ATC"),
    ("1-1", "ATC", "2-1", "ATC"),
    ("2-1", "ATC", "2-1", "RFP"),
    ("2-1", "VCM", "2-1", "RFP"),
    ("2-1", "ATC", "2-1", "VGM"),
]


def main():
    by_replica = pd.read_csv(TABLES / "adsorption_by_replica.csv")
    per_molecule = pd.read_csv(TABLES / "adsorption_per_molecule.csv")
    rows = []

    pooled = by_replica.groupby(["system", "replicate"]).contact_pct_final_window.mean()
    groups = [pooled.loc[s].values for s in SYSTEMS]
    f_stat, p_anova = stats.f_oneway(*groups)
    rows.append({"comparison": "formulation (replica means)", "test": "one-way ANOVA",
                 "n_per_group": 3, "statistic": round(f_stat, 3),
                 "p_value": round(p_anova, 5)})

    molecule_groups = [per_molecule[per_molecule.system == s].contact_pct_final_window
                       for s in SYSTEMS]
    h_stat, p_kw = stats.kruskal(*molecule_groups)
    rows.append({"comparison": "formulation (individual molecules)",
                 "test": "Kruskal-Wallis", "n_per_group": len(molecule_groups[0]),
                 "statistic": round(h_stat, 3), "p_value": f"{p_kw:.3e}"})

    for s1, a1, s2, a2 in PAIRS:
        x = per_molecule[(per_molecule.system == s1) & (per_molecule.antibiotic == a1)]
        y = per_molecule[(per_molecule.system == s2) & (per_molecule.antibiotic == a2)]
        u_stat, p_mw = stats.mannwhitneyu(x.contact_pct_final_window,
                                          y.contact_pct_final_window)
        rows.append({
            "comparison": f"{SYSTEM_LABEL[s1]}/{a1} vs {SYSTEM_LABEL[s2]}/{a2}",
            "test": "Mann-Whitney U", "n_per_group": len(x),
            "statistic": round(u_stat, 1), "p_value": round(p_mw, 4)})

    tests = pd.DataFrame(rows)
    tests.to_csv(TABLES / "statistical_tests.csv", index=False)
    print(tests.to_string(index=False))


if __name__ == "__main__":
    main()
