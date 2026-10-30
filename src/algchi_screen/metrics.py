"""Adsorption metrics computed from direct interatomic contact.

Contact is the primary observable in this study: a molecule counts as adsorbed
in a frame when at least one of its heavy atoms lies within
``CONTACT_CUTOFF_A`` of a polymer heavy atom. Unlike a criterion referenced to
the nanoparticle radius, this is unaffected by particle swelling or dispersion.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .config import (COPIES_PER_ANTIBIOTIC, FINAL_WINDOW_NS, LONG_EVENT_NS,
                     REPLICATES, STABLE_BINDING_FRACTION, TRAJECTORY_NS)


def _final_window(positions):
    return positions["Time_ns"] >= (TRAJECTORY_NS - FINAL_WINDOW_NS)


def contact_fraction(positions, replicate):
    """Fraction of frames in direct contact, whole run and final window."""
    contact = positions[f"Rep{replicate}_Direct_contact"].fillna(0)
    tail = _final_window(positions)
    return float(contact.mean() * 100), float(contact[tail].mean() * 100)


def per_molecule_contact(positions, replicate):
    """Contact fraction of each individual molecule over the final window."""
    tail = _final_window(positions)
    grouped = positions[tail].groupby("Molecule_ID")[f"Rep{replicate}_Direct_contact"]
    return grouped.mean() * 100


def uptake(positions, replicate):
    """Molecules that ever contacted, and molecules stably bound at the end."""
    col = f"Rep{replicate}_Direct_contact"
    ever = int((positions.groupby("Molecule_ID")[col].max() == 1).sum())
    tail = _final_window(positions)
    stable = int((positions[tail].groupby("Molecule_ID")[col].mean()
                  > STABLE_BINDING_FRACTION).sum())
    return ever, stable


def _runs(binary):
    """Lengths of maximal runs of consecutive 1s."""
    out, count = [], 0
    for value in binary:
        if value == 1:
            count += 1
        else:
            if count:
                out.append(count)
            count = 0
    if count:
        out.append(count)
    return out


def contact_events(positions, replicate):
    """Durations in ns of every contact event, pooled over the eight molecules."""
    col = f"Rep{replicate}_Direct_contact"
    durations = []
    for _, group in positions.groupby("Molecule_ID"):
        group = group.sort_values("Frame")
        step = group["Time_ns"].diff().median()
        durations.extend(length * step for length in
                         _runs(group[col].fillna(0).astype(int).values))
    return np.array(durations)


def summarise_condition(positions, system, antibiotic):
    """All contact metrics for one system-antibiotic, one row per replica."""
    rows = []
    for replicate in REPLICATES:
        whole, final = contact_fraction(positions, replicate)
        ever, stable = uptake(positions, replicate)
        events = contact_events(positions, replicate)
        rows.append({
            "system": system,
            "antibiotic": antibiotic,
            "replicate": replicate,
            "n_molecules": COPIES_PER_ANTIBIOTIC,
            "contact_pct_whole_run": whole,
            "contact_pct_final_window": final,
            "molecules_ever_contacting": ever,
            "molecules_stably_bound": stable,
            "n_contact_events": len(events),
            "median_residence_ns": float(np.median(events)) if len(events) else 0.0,
            "max_residence_ns": float(events.max()) if len(events) else 0.0,
            "events_over_10ns": int((events >= LONG_EVENT_NS).sum()),
        })
    return pd.DataFrame(rows)


def classifier_states(state_counts, replicate):
    """Mean occupancy of the radius-referenced FREE/ADSORBED/ABSORBED states.

    Reported only for comparison with the contact criterion; see the manuscript
    for why this definition fails on a particle whose boundary moves.
    """
    free = state_counts[f"Rep{replicate}_Free"]
    adsorbed = state_counts[f"Rep{replicate}_Adsorbed"]
    absorbed = state_counts[f"Rep{replicate}_Absorbed"]
    total = free + adsorbed + absorbed
    return {
        "free_pct": float((free / total * 100).mean()),
        "adsorbed_pct": float((adsorbed / total * 100).mean()),
        "absorbed_pct": float((absorbed / total * 100).mean()),
        "bound_pct": float(((adsorbed + absorbed) / total * 100).mean()),
    }


def absorbed_frame_counts(positions):
    """Frames classified ABSORBED, and how many of those had no polymer contact."""
    total = without_contact = 0
    for replicate in REPLICATES:
        absorbed = positions[f"Rep{replicate}_State"] == "ABSORBED"
        total += int(absorbed.sum())
        without_contact += int(
            (absorbed & (positions[f"Rep{replicate}_Direct_contact"] == 0)).sum())
    return total, without_contact


def absorbed_without_contact(positions):
    """Percentage of ABSORBED frames in which no polymer atom is within the cutoff.

    A large value means the classification boundary expanded past molecules that
    never touched the polymer, i.e. geometric engulfment rather than uptake.
    Frames are pooled across replicas, so the statistic weights every classified
    frame equally.
    """
    total, without_contact = absorbed_frame_counts(positions)
    return 100.0 * without_contact / total if total else 0.0
