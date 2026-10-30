"""Separation of polymer-internal from antibiotic-specific interactions.

Every antibiotic in a given system shares one trajectory, so the interaction
totals recorded against it contain a large contribution from the polymer's own
internal network. Taking the frame-wise minimum across the antibiotics of a
system estimates that shared contribution: subtracting it leaves the excess
attributable to a particular antibiotic, while the minimum itself describes the
polymer network.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from .config import CHAIN_COUNTS, PANEL, REPLICATES
from .io import load_interactions

COMPONENTS = ("Alginate", "Chitosan", "TotalPolymer")


def _stack(system, component, replicate, metric):
    frames = [load_interactions(system, ab, component)[f"Rep{replicate}_{metric}"]
              .values.astype(float) for ab in PANEL[system]]
    return np.vstack(frames)


def _nan_safe(function, array, axis=None):
    """Apply a nan-aware reduction, tolerating columns that are empty by design.

    Some interaction types (for example pi-cation with alginate) are absent for
    a whole trajectory and are stored as all-NaN; numpy warns on those slices
    even though an empty result is the correct answer.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = function(array, axis=axis) if axis is not None else function(array)
    return np.nan_to_num(result)


def polymer_baseline(system, component, replicate, metric="Hydrogen_bonds"):
    """Frame-wise minimum across antibiotics: the polymer-internal contribution."""
    return _nan_safe(np.nanmin, _stack(system, component, replicate, metric), axis=0)


def antibiotic_excess(system, antibiotic, component, replicate,
                      metric="Hydrogen_bonds"):
    """Counts in excess of the polymer-internal baseline for one antibiotic."""
    stack = _stack(system, component, replicate, metric)
    index = PANEL[system].index(antibiotic)
    return np.nan_to_num(stack[index]) - _nan_safe(np.nanmin, stack, axis=0)


def network_decomposition(system):
    """Polymer network per replica, split into its three interaction classes.

    The alginate-chitosan interfacial term is obtained by difference, since the
    total-polymer count includes it while the component counts do not.
    """
    counts = CHAIN_COUNTS[system]
    rows = []
    for replicate in REPLICATES:
        alginate = polymer_baseline(system, "Alginate", replicate)
        chitosan = polymer_baseline(system, "Chitosan", replicate)
        total = polymer_baseline(system, "TotalPolymer", replicate)
        interfacial = total - alginate - chitosan
        salt = polymer_baseline(system, "TotalPolymer", replicate, "Salt_bridges")
        rows.append({
            "system": system,
            "replicate": replicate,
            "alginate_alginate_hbonds": float(_nan_safe(np.nanmean, alginate)),
            "chitosan_chitosan_hbonds": float(_nan_safe(np.nanmean, chitosan)),
            "alginate_chitosan_hbonds": float(_nan_safe(np.nanmean, interfacial)),
            "total_hbonds": float(_nan_safe(np.nanmean, total)),
            "hbonds_per_alginate_chain": float(_nan_safe(np.nanmean, alginate)) / counts["alginate"],
            "interfacial_hbonds_per_chitosan_chain":
                float(_nan_safe(np.nanmean, interfacial)) / counts["chitosan"],
            "calcium_salt_bridges": float(_nan_safe(np.nanmean, salt)),
            "salt_bridges_per_alginate_chain":
                float(_nan_safe(np.nanmean, salt)) / counts["alginate"],
        })
    return pd.DataFrame(rows)


def antibiotic_excess_summary(system, antibiotic):
    """Antibiotic-specific hydrogen bonds per frame, one row per replica."""
    rows = []
    for replicate in REPLICATES:
        excess = antibiotic_excess(system, antibiotic, "TotalPolymer", replicate)
        chitosan = antibiotic_excess(system, antibiotic, "Chitosan", replicate)
        pication = antibiotic_excess(system, antibiotic, "TotalPolymer",
                                     replicate, "Pi_cation")
        rows.append({
            "system": system,
            "antibiotic": antibiotic,
            "replicate": replicate,
            "excess_hbonds_total_polymer": float(_nan_safe(np.nanmean, excess)),
            "excess_hbonds_chitosan": float(_nan_safe(np.nanmean, chitosan)),
            "excess_pi_cation": float(_nan_safe(np.nanmean, pication)),
        })
    return pd.DataFrame(rows)
