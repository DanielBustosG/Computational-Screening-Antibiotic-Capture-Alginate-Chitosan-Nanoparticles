"""02 - Nanoparticle structural stability.

Inputs : data/raw/md/NP_geometry_consolidado.xlsx, Rg_consolidado.xlsx,
         SASA_consolidado.xlsx, Interactions_*.xlsx
Outputs: results/tables/stability_by_replica.csv
         results/tables/stability_summary.csv   (manuscript Table S4)
         results/tables/polymer_network.csv

Reproduces the finding that the 2:1 particle is structurally invariant, the 1:1
particle swells moderately, and the 1:2 particle disperses.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from algchi_screen.config import CHAIN_COUNTS, SYSTEMS, SYSTEM_LABEL, TABLES
from algchi_screen.interactions import network_decomposition
from algchi_screen.io import load_geometry, load_rg, load_sasa
from algchi_screen.stability import stability_by_replica


def main():
    TABLES.mkdir(parents=True, exist_ok=True)

    per_replica = pd.concat(
        [stability_by_replica(load_geometry(s), load_rg(s), load_sasa(s), s)
         for s in SYSTEMS],
        ignore_index=True)
    per_replica.to_csv(TABLES / "stability_by_replica.csv", index=False)

    network = pd.concat([network_decomposition(s) for s in SYSTEMS],
                        ignore_index=True)
    network.to_csv(TABLES / "polymer_network.csv", index=False)

    rows = []
    for system in SYSTEMS:
        geo = per_replica[per_replica.system == system]
        net = network[network.system == system]
        counts = CHAIN_COUNTS[system]
        rows.append({
            "Formulation": SYSTEM_LABEL[system],
            "Alginate chains": counts["alginate"],
            "Chitosan chains": counts["chitosan"],
            "Calcium ions": counts["calcium"],
            "Initial radius P95 (A)": round(geo.initial_radius_A.mean(), 1),
            "Final radius P95 (A)": f"{geo.final_radius_A.mean():.1f} +/- {geo.final_radius_A.std():.1f}",
            "Radial compaction": f"{geo.radial_compaction.mean():.2f} +/- {geo.radial_compaction.std():.2f}",
            "Apparent volume ratio": f"{geo.apparent_volume_ratio.mean():.2f} +/- {geo.apparent_volume_ratio.std():.2f}",
            "Final Rg (A)": f"{geo.final_Rg_A.mean():.1f} +/- {geo.final_Rg_A.std():.1f}",
            "Initial SASA (1e3 A2)": round(geo.initial_SASA_A2.mean() / 1000),
            "Final SASA (1e3 A2)": f"{geo.final_SASA_A2.mean()/1000:.0f} +/- {geo.final_SASA_A2.std()/1000:.0f}",
            "Shape anisotropy": round(geo.mean_shape_anisotropy.mean(), 3),
            "Alg-Alg H-bonds per alginate chain": f"{net.hbonds_per_alginate_chain.mean():.2f} +/- {net.hbonds_per_alginate_chain.std():.2f}",
            "Alg-Chi H-bonds per chitosan chain": f"{net.interfacial_hbonds_per_chitosan_chain.mean():.2f} +/- {net.interfacial_hbonds_per_chitosan_chain.std():.2f}",
            "Calcium salt bridges per frame": f"{net.calcium_salt_bridges.mean():.0f} +/- {net.calcium_salt_bridges.std():.0f}",
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(TABLES / "stability_summary.csv", index=False)
    summary.to_excel(TABLES / "stability_summary.xlsx", index=False)

    print("Apparent volume ratio after 100 ns (manuscript Section 3.7):")
    for _, row in summary.iterrows():
        print(f"  {row['Formulation']}: {row['Apparent volume ratio']}")
    print("High-asphericity frames across all systems and replicas: "
          f"{int(per_replica.frames_high_asphericity.sum())}")


if __name__ == "__main__":
    main()
