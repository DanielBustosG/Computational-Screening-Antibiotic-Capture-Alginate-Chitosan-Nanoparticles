"""05 - Manuscript figures 5 and 6 and supplementary figures S4 and S5.

Inputs : data/raw/md/*.xlsx, results/tables/adsorption_*.csv
Outputs: results/figures/Fig5.png   (nanoparticle stability)
         results/figures/Fig6.png   (antibiotic adsorption)
         results/figures/FigS4.png  (density profiles)
         results/figures/FigS5.png  (radius-classifier artifact)
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from algchi_screen.config import (ANTIBIOTIC_COLOR, FIGURE_DPI, FIGURES, PANEL,
                                  REPLICATES, SYSTEM_COLOR, SYSTEM_LABEL,
                                  SYSTEMS, TABLES)
from algchi_screen.interactions import network_decomposition, polymer_baseline
from algchi_screen.io import (load_density, load_diffusion, load_geometry,
                              load_positions, load_rg)
from algchi_screen.metrics import absorbed_frame_counts

ORDER = [(s, a) for s in SYSTEMS for a in PANEL[s]]

plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9,
                     "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 110})


def _band(ax, time, values, color):
    mean = np.nanmean(values, axis=0)
    sd = np.nanstd(values, axis=0)
    ax.plot(time, mean, color=color, lw=1.4)
    ax.fill_between(time, mean - sd, mean + sd, color=color, alpha=0.2, lw=0)
    return mean


def figure_stability():
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.6))
    for system in SYSTEMS:
        geo = load_geometry(system)
        stack = np.vstack([geo[f"Rep{r}_NP_radius_P95_A"].values.astype(float)
                           for r in REPLICATES])
        _band(axes[0], geo["Time_ns"], stack, SYSTEM_COLOR[system])
        rg = load_rg(system)
        stack = np.vstack([rg[f"Rep{r}_Rg_A"].values.astype(float) for r in REPLICATES])
        _band(axes[1], rg["Time"], stack, SYSTEM_COLOR[system])

    axes[0].set(xlabel="Time (ns)", ylabel="Particle radius, P95 (Å)")
    axes[0].set_title("Chitosan-rich particle disperses", pad=6)
    axes[1].set(xlabel="Time (ns)", ylabel="Radius of gyration (Å)")
    axes[1].set_title("Expansion is a real swelling", pad=6)
    for ax, heights in ((axes[0], (131, 67, 52)), (axes[1], (62, 41.4, 38.6))):
        for label, system, y in zip(("1:2", "1:1", "2:1"), ("1-2", "1-1", "2-1"), heights):
            ax.text(101, y, label, color=SYSTEM_COLOR[system], fontsize=8,
                    fontweight="bold", va="center")
        ax.set_xlim(-3, 118)

    network = pd.concat([network_decomposition(s) for s in SYSTEMS], ignore_index=True)
    x = np.arange(3)
    alg = [network[network.system == s].hbonds_per_alginate_chain for s in SYSTEMS]
    interfacial = [network[network.system == s].interfacial_hbonds_per_chitosan_chain
                   for s in SYSTEMS]
    axes[2].bar(x - 0.19, [v.mean() for v in alg], 0.36,
                yerr=[v.std() for v in alg], color="#4b6584", capsize=2,
                label="alginate-alginate\nper alginate chain")
    axes[2].bar(x + 0.19, [v.mean() for v in interfacial], 0.36,
                yerr=[v.std() for v in interfacial], color="#c9a227", capsize=2,
                label="alginate-chitosan\nper chitosan chain")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels([SYSTEM_LABEL[s] for s in SYSTEMS])
    axes[2].set(xlabel="Alginate:chitosan ratio", ylabel="Hydrogen bonds per chain",
                ylim=(0, 8.2))
    axes[2].set_title("Cross-link density explains it", pad=6)
    axes[2].legend(frameon=False, fontsize=6.5, loc="upper left",
                   handlelength=1.1, borderpad=0.1)
    for letter, ax in zip("abc", axes):
        ax.text(-0.20, 1.10, letter, transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top")
    fig.tight_layout(w_pad=1.8)
    fig.savefig(FIGURES / "Fig5.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def figure_adsorption():
    per_molecule = pd.read_csv(TABLES / "adsorption_per_molecule.csv")
    by_replica = pd.read_csv(TABLES / "adsorption_by_replica.csv")
    events = pd.read_csv(TABLES / "contact_events.csv")
    diffusion = load_diffusion()
    diffusion["D_mean"] = diffusion[[f"Rep{r}_D_fit_A2_per_ns" for r in REPLICATES]].mean(axis=1)
    diffusion["D_sd"] = diffusion[[f"Rep{r}_D_fit_A2_per_ns" for r in REPLICATES]].std(axis=1)

    fig = plt.figure(figsize=(7.4, 5.3))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.92], hspace=0.60, wspace=0.44)
    bands = [(-0.5, 2.5, "1-1"), (2.5, 5.5, "1-2"), (5.5, 9.5, "2-1")]

    ax = fig.add_subplot(gs[0, :2])
    for x0, x1, system in bands:
        ax.axvspan(x0, x1, color=SYSTEM_COLOR[system], alpha=0.06, lw=0)
        ax.text((x0 + x1) / 2, 106, SYSTEM_LABEL[system], ha="center", fontsize=8,
                color=SYSTEM_COLOR[system], fontweight="bold")
    for i, (system, antibiotic) in enumerate(ORDER):
        values = per_molecule[(per_molecule.system == system)
                              & (per_molecule.antibiotic == antibiotic)
                              ].contact_pct_final_window.values
        jitter = np.random.default_rng(i).normal(i, 0.10, len(values))
        ax.scatter(jitter, values, s=7, color=ANTIBIOTIC_COLOR[antibiotic],
                   alpha=0.55, lw=0)
        ax.plot([i - 0.30, i + 0.30], [values.mean()] * 2,
                color=ANTIBIOTIC_COLOR[antibiotic], lw=2.2, solid_capstyle="butt")
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([a for _, a in ORDER], fontsize=7)
    ax.set(ylabel="Time in direct contact (%)", ylim=(-6, 116), xlim=(-0.6, 9.6))
    ax.set_title("Adsorption from direct polymer contact", pad=8)

    ax = fig.add_subplot(gs[0, 2])
    for i, system in enumerate(SYSTEMS):
        subset = by_replica[by_replica.system == system]
        breadth = subset.groupby("replicate").molecules_ever_contacting.mean() / 8 * 100
        ax.bar(i, breadth.mean(), 0.62, yerr=breadth.std(),
               color=SYSTEM_COLOR[system], capsize=2.5)
        ax.scatter([i] * len(breadth), breadth.values, s=9, color="0.15", zorder=3)
    ax.set_xticks(range(3))
    ax.set_xticklabels([SYSTEM_LABEL[s] for s in SYSTEMS])
    ax.set(ylabel="Molecules contacted (%)", xlabel="Alg:Chi ratio", ylim=(0, 122))
    ax.set_title("Uptake breadth", pad=8)

    ax = fig.add_subplot(gs[1, 0])
    # The three survival curves nearly superimpose, so a legend is used instead of
    # inline labels, which would overlap and misattribute curves.
    for system in SYSTEMS:
        durations = np.sort(events[events.system == system].duration_ns.values)
        survival = 1 - np.arange(len(durations)) / len(durations)
        ax.step(durations, survival, color=SYSTEM_COLOR[system], lw=1.5, where="post",
                label=f"{SYSTEM_LABEL[system]}  (n = {len(durations)})")
    ax.legend(frameon=False, fontsize=6.3, loc="lower left", handlelength=1.2,
              borderpad=0.1, labelspacing=0.25)
    ax.set(xscale="log", yscale="log", xlabel="Contact duration (ns)",
           ylabel="P(duration > t)", xlim=(0.35, 140), ylim=(3e-3, 1.9))
    ax.set_title("Residence times", pad=8)

    ax = fig.add_subplot(gs[1, 1])
    for x0, x1, system in bands:
        ax.axvspan(x0, x1, color=SYSTEM_COLOR[system], alpha=0.06, lw=0)
    for i, (system, antibiotic) in enumerate(ORDER):
        row = diffusion[(diffusion.system == system)
                        & (diffusion.antibiotic == antibiotic)]
        ax.errorbar(i, row.D_mean.iloc[0], yerr=row.D_sd.iloc[0], fmt="o", ms=4,
                    color=ANTIBIOTIC_COLOR[antibiotic], capsize=2, lw=1.1)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([a for _, a in ORDER], rotation=90, fontsize=6.5)
    ax.set(ylabel="D (Å² ns⁻¹)", xlim=(-0.6, 9.6))
    ax.set_title("Diffusion", pad=8)

    ax = fig.add_subplot(gs[1, 2])
    for x0, x1, system in bands:
        ax.axvspan(x0, x1, color=SYSTEM_COLOR[system], alpha=0.06, lw=0)
    summary = pd.read_csv(TABLES / "adsorption_by_replica.csv")
    from algchi_screen.interactions import antibiotic_excess_summary
    for i, (system, antibiotic) in enumerate(ORDER):
        excess = antibiotic_excess_summary(system, antibiotic)
        ax.bar(i, excess.excess_hbonds_total_polymer.mean(), 0.72,
               yerr=excess.excess_hbonds_total_polymer.std(),
               color=ANTIBIOTIC_COLOR[antibiotic], capsize=1.8,
               error_kw=dict(lw=0.9))
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([a for _, a in ORDER], rotation=90, fontsize=6.5)
    ax.set(ylabel="Excess H-bonds (per frame)", xlim=(-0.6, 9.6))
    ax.set_title("Antibiotic-specific bonding", pad=8)

    for letter, ax in zip("abcde", fig.axes):
        ax.text(-0.20, 1.12, letter, transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top")
    fig.savefig(FIGURES / "Fig6.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def figure_classifier_artifact():
    comparison = pd.read_csv(TABLES / "classifier_vs_contact.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.6))
    x = np.arange(3)
    labels = [SYSTEM_LABEL[s] for s in SYSTEMS]
    radius = [comparison[comparison.Formulation == l]["Radius-based bound (%)"].mean()
              for l in labels]
    radius_sd = [comparison[comparison.Formulation == l]["Radius-based bound (%)"].std()
                 for l in labels]
    contact = [comparison[comparison.Formulation == l]["Contact-based bound (%)"].mean()
               for l in labels]
    contact_sd = [comparison[comparison.Formulation == l]["Contact-based bound (%)"].std()
                  for l in labels]
    axes[0].bar(x - 0.175, radius, 0.35, yerr=radius_sd, color="#b0b0b0", capsize=2,
                label="radius-based state\n(ADSORBED + ABSORBED)")
    axes[0].bar(x + 0.175, contact, 0.35, yerr=contact_sd, color="#2b3a67", capsize=2,
                label="direct contact ≤ 3.5 Å")
    axes[0].set_xticks(x); axes[0].set_xticklabels(labels)
    axes[0].set(xlabel="Alg:Chi ratio", ylabel="Antibiotic bound (%)", ylim=(0, 118))
    axes[0].set_title("Two definitions disagree in 1:2", pad=6)
    axes[0].legend(frameon=False, fontsize=6.3, loc="upper left", handlelength=1.1)

    fractions = []
    for system in SYSTEMS:
        total = without = 0
        for antibiotic in PANEL[system]:
            t, w = absorbed_frame_counts(load_positions(system, antibiotic))
            total += t
            without += w
        fractions.append(100.0 * without / total if total else 0.0)
    axes[1].bar(x, fractions, 0.6, color=[SYSTEM_COLOR[s] for s in SYSTEMS])
    for i, value in enumerate(fractions):
        axes[1].text(i, value + 2.5, f"{value:.0f}%", ha="center", fontsize=7)
    axes[1].set_xticks(x); axes[1].set_xticklabels(labels)
    axes[1].set(xlabel="Alg:Chi ratio",
                ylabel="ABSORBED frames without\nany polymer contact (%)", ylim=(0, 80))
    axes[1].set_title("Engulfment, not absorption", pad=6)

    geo = load_geometry("1-2")
    positions = load_positions("1-2", "ATC")
    axes[2].plot(geo["Time_ns"], geo["Rep1_Core_radius_A"], color="#c8622a", lw=1.5,
                 label="core radius (classifier)")
    axes[2].plot(geo["Time_ns"], geo["Rep1_NP_radius_P95_A"], color="#c8622a", lw=1.0,
                 ls="--", label="outer radius P95")
    median = positions.groupby("Time_ns")["Rep1_COM_distance_A"].median()
    axes[2].plot(median.index, median.values, color="#2b3a67", lw=1.3,
                 label="median antibiotic\ndistance from centre")
    axes[2].set(xlabel="Time (ns)", ylabel="Distance from centre (Å)", ylim=(0, 175))
    axes[2].set_title("Boundary overtakes the ligands", pad=6)
    axes[2].legend(frameon=False, fontsize=6.3, loc="lower right", handlelength=1.3)
    for letter, ax in zip("abc", axes):
        ax.text(-0.20, 1.10, letter, transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top")
    fig.tight_layout(w_pad=1.6)
    fig.savefig(FIGURES / "FigS5.png",
                dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def figure_density():
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.5), sharex=True, sharey=True)
    for ax, system in zip(axes, SYSTEMS):
        polymer = load_density(system, "Polymer")["Z"]
        antibiotics = load_density(system, "Antibiotics")["Z"]
        x = polymer["Position_A"].values.astype(float)
        alginate = np.nanmean([np.nan_to_num(polymer[f"ALG_rep{r}_mean"].values.astype(float))
                               for r in REPLICATES], axis=0)
        centre = np.trapezoid(alginate * x, x) / np.trapezoid(alginate, x)
        for tag, color, label in (("ALG", "#2b3a67", "alginate"),
                                  ("CHI", "#c9a227", "chitosan")):
            profile = np.nanmean([np.nan_to_num(polymer[f"{tag}_rep{r}_mean"].values.astype(float))
                                  for r in REPLICATES], axis=0)
            profile = profile / np.nanmax(profile)
            ax.fill_between(x - centre, 0, profile, color=color, alpha=0.35, lw=0)
            ax.plot(x - centre, profile, color=color, lw=1.1, label=label)
        xa = antibiotics["Position_A"].values.astype(float)
        combined = np.nanmean(
            [np.nan_to_num(antibiotics[f"{a}_rep{r}_mean"].values.astype(float))
             for a in PANEL[system] for r in REPLICATES], axis=0)
        ax.plot(xa - centre, combined / np.nanmax(combined), color="#d62728", lw=1.3,
                label="antibiotics")
        ax.set_title(SYSTEM_LABEL[system], pad=6)
    axes[0].set_ylabel("Normalised density")
    axes[0].legend(frameon=False, fontsize=6.5, loc="upper left",
                   handlelength=1.2, borderpad=0.1)
    axes[1].set_xlabel("Distance from particle centre, z (Å)")
    axes[-1].set(xlim=(-150, 150), ylim=(0, 1.25))
    for letter, ax in zip("abc", axes):
        ax.text(-0.16, 1.12, letter, transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top")
    fig.tight_layout(w_pad=0.8)
    fig.savefig(FIGURES / "FigS4.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    figure_stability()
    figure_adsorption()
    figure_classifier_artifact()
    figure_density()
    for path in sorted(FIGURES.glob("*.png")):
        print(f"{path.name}: {path.stat().st_size / 1024:.0f} kB")


if __name__ == "__main__":
    main()
