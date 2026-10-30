"""03 - Antibiotic adsorption from direct interatomic contact.

Inputs : data/raw/md/Classifier_positions_consolidado.xlsx,
         Classifier_state_counts_consolidado.xlsx, MSD_DIFF_consolidado.xlsx,
         Interactions_*.xlsx
Outputs: results/tables/adsorption_by_replica.csv
         results/tables/adsorption_per_molecule.csv
         results/tables/adsorption_summary.csv     (manuscript Table S3)
         results/tables/contact_events.csv
         results/tables/classifier_vs_contact.csv  (manuscript Table S5)

Adsorption is measured by direct contact rather than by radial position, so the
numbers are unaffected by the expansion of the 1:2 particle.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from algchi_screen.config import (ANTIBIOTIC_NAME, PANEL, REPLICATES, SYSTEMS,
                                  SYSTEM_LABEL, TABLES)
from algchi_screen.interactions import antibiotic_excess_summary
from algchi_screen.io import load_diffusion, load_positions, load_state_counts
from algchi_screen.metrics import (absorbed_frame_counts, absorbed_without_contact,
                                   classifier_states, contact_events,
                                   per_molecule_contact, summarise_condition)


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    by_replica, per_molecule, events, excess, comparison = [], [], [], [], []

    for system in SYSTEMS:
        for antibiotic in PANEL[system]:
            positions = load_positions(system, antibiotic)
            by_replica.append(summarise_condition(positions, system, antibiotic))
            excess.append(antibiotic_excess_summary(system, antibiotic))

            for replicate in REPLICATES:
                molecules = per_molecule_contact(positions, replicate)
                per_molecule.append(pd.DataFrame({
                    "system": system, "antibiotic": antibiotic,
                    "replicate": replicate,
                    "molecule": molecules.index,
                    "contact_pct_final_window": molecules.values}))
                durations = contact_events(positions, replicate)
                events.append(pd.DataFrame({
                    "system": system, "antibiotic": antibiotic,
                    "replicate": replicate, "duration_ns": durations}))

            states = load_state_counts(system, antibiotic)
            state_means = pd.DataFrame(
                [classifier_states(states, r) for r in REPLICATES]).mean()
            contact_mean = by_replica[-1].contact_pct_whole_run.mean()
            comparison.append({
                "Formulation": SYSTEM_LABEL[system],
                "Antibiotic": ANTIBIOTIC_NAME[antibiotic],
                "FREE (%)": round(state_means.free_pct, 1),
                "ADSORBED (%)": round(state_means.adsorbed_pct, 1),
                "ABSORBED (%)": round(state_means.absorbed_pct, 1),
                "Radius-based bound (%)": round(state_means.bound_pct, 1),
                "Contact-based bound (%)": round(contact_mean, 1),
                "Discrepancy (pp)": round(state_means.bound_pct - contact_mean, 1),
                "ABSORBED frames with no contact (%)":
                    round(absorbed_without_contact(positions), 1),
            })

    by_replica = pd.concat(by_replica, ignore_index=True)
    per_molecule = pd.concat(per_molecule, ignore_index=True)
    events = pd.concat(events, ignore_index=True)
    excess = pd.concat(excess, ignore_index=True)
    comparison = pd.DataFrame(comparison)

    by_replica.to_csv(TABLES / "adsorption_by_replica.csv", index=False)
    per_molecule.to_csv(TABLES / "adsorption_per_molecule.csv", index=False)
    events.to_csv(TABLES / "contact_events.csv", index=False)
    comparison.to_csv(TABLES / "classifier_vs_contact.csv", index=False)
    comparison.to_excel(TABLES / "classifier_vs_contact.xlsx", index=False)

    diffusion = load_diffusion()
    diffusion["D_mean_A2_per_ns"] = diffusion[
        [f"Rep{r}_D_fit_A2_per_ns" for r in REPLICATES]].mean(axis=1)
    diffusion["D_sd_A2_per_ns"] = diffusion[
        [f"Rep{r}_D_fit_A2_per_ns" for r in REPLICATES]].std(axis=1)

    rows = []
    for system in SYSTEMS:
        for antibiotic in PANEL[system]:
            reps = by_replica[(by_replica.system == system)
                              & (by_replica.antibiotic == antibiotic)]
            exc = excess[(excess.system == system)
                         & (excess.antibiotic == antibiotic)]
            dif = diffusion[(diffusion.system == system)
                            & (diffusion.antibiotic == antibiotic)]
            rows.append({
                "Formulation": SYSTEM_LABEL[system],
                "Antibiotic": ANTIBIOTIC_NAME[antibiotic],
                "Molecules contacted (of 8)":
                    f"{reps.molecules_ever_contacting.mean():.1f} +/- {reps.molecules_ever_contacting.std():.1f}",
                "Contact time, whole run (%)":
                    f"{reps.contact_pct_whole_run.mean():.1f} +/- {reps.contact_pct_whole_run.std():.1f}",
                "Contact time, final 20 ns (%)":
                    f"{reps.contact_pct_final_window.mean():.1f} +/- {reps.contact_pct_final_window.std():.1f}",
                "Stably bound molecules (of 8)":
                    f"{reps.molecules_stably_bound.mean():.1f} +/- {reps.molecules_stably_bound.std():.1f}",
                "Contact events (per run)":
                    f"{reps.n_contact_events.mean():.0f} +/- {reps.n_contact_events.std():.0f}",
                "Median residence (ns)": round(reps.median_residence_ns.mean(), 2),
                "Longest contact (ns)":
                    f"{reps.max_residence_ns.mean():.1f} +/- {reps.max_residence_ns.std():.1f}",
                "Events > 10 ns": round(reps.events_over_10ns.mean(), 1),
                "Excess H-bonds (per frame)":
                    f"{exc.excess_hbonds_total_polymer.mean():.2f} +/- {exc.excess_hbonds_total_polymer.std():.2f}",
                "D (A2/ns)":
                    f"{dif.D_mean_A2_per_ns.iloc[0]:.1f} +/- {dif.D_sd_A2_per_ns.iloc[0]:.1f}",
            })
    summary = pd.DataFrame(rows)
    summary.to_csv(TABLES / "adsorption_summary.csv", index=False)
    summary.to_excel(TABLES / "adsorption_summary.xlsx", index=False)

    print("Contact time over the final 20 ns, pooled by formulation "
          "(manuscript Section 3.8):")
    pooled = by_replica.groupby(["system", "replicate"]).contact_pct_final_window.mean()
    for system in SYSTEMS:
        values = pooled.loc[system]
        print(f"  {SYSTEM_LABEL[system]}: {values.mean():.1f} +/- {values.std():.1f} %")
    print("ABSORBED frames without any polymer contact, pooled per formulation:")
    for system in SYSTEMS:
        total = without = 0
        for antibiotic in PANEL[system]:
            t, w = absorbed_frame_counts(load_positions(system, antibiotic))
            total += t
            without += w
        share = 100.0 * without / total if total else 0.0
        print(f"  {SYSTEM_LABEL[system]}: {share:.1f} % ({without} of {total} frames)")


if __name__ == "__main__":
    main()
