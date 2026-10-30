"""Loading and reshaping of the consolidated trajectory workbooks.

The workbooks store one sheet per system-antibiotic (or per system) with the
three replicas side by side as ``Rep1_``/``Rep2_``/``Rep3_`` column prefixes.
These helpers convert that wide layout into long dataframes.
"""
from __future__ import annotations

import pandas as pd

from .config import MD_DIR, REPLICATES


def _sheet(path, name):
    return pd.read_excel(path, sheet_name=name)


def load_positions(system, antibiotic):
    """Per-molecule, per-frame positional and contact data for one condition."""
    path = MD_DIR / "Classifier_positions_consolidado.xlsx"
    return _sheet(path, f"{system}_{antibiotic}")


def load_state_counts(system, antibiotic):
    """Per-frame counts of free / adsorbed / absorbed molecules."""
    path = MD_DIR / "Classifier_state_counts_consolidado.xlsx"
    return _sheet(path, f"{system}_{antibiotic}")


def load_geometry(system):
    """Per-frame nanoparticle radius, radius of gyration and shape descriptors."""
    path = MD_DIR / "NP_geometry_consolidado.xlsx"
    return _sheet(path, f"System_{system}")


def load_rg(system):
    path = MD_DIR / "Rg_consolidado.xlsx"
    return _sheet(path, f"System_{system}")


def load_sasa(system):
    path = MD_DIR / "SASA_consolidado.xlsx"
    return _sheet(path, f"System_{system}")


def load_interactions(system, antibiotic, component):
    """Interaction counts per frame.

    ``component`` is one of ``Alginate``, ``Chitosan`` or ``TotalPolymer``.
    These totals include the polymer's own internal interactions; use
    :func:`algchi_screen.interactions.polymer_baseline` or
    :func:`algchi_screen.interactions.antibiotic_excess` to separate them.
    """
    path = MD_DIR / f"Interactions_{system}.xlsx"
    return _sheet(path, f"{antibiotic}_{component}")


def load_diffusion():
    """Fitted diffusion coefficients for every system-antibiotic combination."""
    path = MD_DIR / "MSD_DIFF_consolidado.xlsx"
    fit = _sheet(path, "Fit_summary")
    fit[["system", "antibiotic"]] = fit["System_Antibiotic"].str.split("_", expand=True)
    return fit


def load_density(system, kind):
    """Density profiles. ``kind`` is ``Polymer`` or ``Antibiotics``."""
    path = MD_DIR / f"Density_{system}.xlsx"
    return {axis: _sheet(path, f"{kind}_{axis}") for axis in ("X", "Y", "Z")}


def replica_columns(frame, suffix):
    """Return the per-replica column names for a metric, e.g. ``NP_Rg_A``."""
    return [f"Rep{r}_{suffix}" for r in REPLICATES if f"Rep{r}_{suffix}" in frame.columns]
