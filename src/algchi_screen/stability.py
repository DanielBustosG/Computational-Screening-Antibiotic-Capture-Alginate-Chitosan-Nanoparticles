"""Nanoparticle structural stability descriptors."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .config import CHAIN_COUNTS, FINAL_WINDOW_NS, REPLICATES, TRAJECTORY_NS


def _tail(frame, time_column):
    return frame[time_column] >= (TRAJECTORY_NS - FINAL_WINDOW_NS)


def stability_by_replica(geometry, rg, sasa, system):
    """Initial and final size descriptors for each replica of one system."""
    rows = []
    for replicate in REPLICATES:
        radius = geometry[f"Rep{replicate}_NP_radius_P95_A"]
        tail_geo = _tail(geometry, "Time_ns")
        tail_rg = _tail(rg, "Time")
        tail_sasa = _tail(sasa, "Time")
        final_radius = float(radius[tail_geo].mean())
        initial_radius = float(radius.iloc[0])
        rows.append({
            "system": system,
            "replicate": replicate,
            "initial_radius_A": initial_radius,
            "final_radius_A": final_radius,
            "radial_compaction": final_radius / initial_radius,
            "apparent_volume_ratio": (final_radius / initial_radius) ** 3,
            "initial_Rg_A": float(rg[f"Rep{replicate}_Rg_A"].iloc[0]),
            "final_Rg_A": float(rg.loc[tail_rg, f"Rep{replicate}_Rg_A"].mean()),
            "initial_SASA_A2": float(sasa[f"Rep{replicate}_SASA_A2"].iloc[:10].mean()),
            "final_SASA_A2": float(sasa.loc[tail_sasa, f"Rep{replicate}_SASA_A2"].mean()),
            "mean_shape_anisotropy": float(geometry[f"Rep{replicate}_Shape_anisotropy"].mean()),
            "frames_high_asphericity": int(
                geometry[f"Rep{replicate}_Shape_warning"].astype(str)
                .str.contains("HIGH").sum()),
        })
    return pd.DataFrame(rows)
