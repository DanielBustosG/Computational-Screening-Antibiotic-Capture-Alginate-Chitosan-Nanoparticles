"""Paths, frozen parameters and input checksums for the analysis.

Every script imports from here; no path or threshold is duplicated elsewhere.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA_RAW = ROOT / "data" / "raw"
DATA_PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"
FIGURES = RESULTS / "figures"

MD_DIR = DATA_RAW / "md"
DOCK_DIR = DATA_RAW / "docking"

# --- Simulated systems -------------------------------------------------------
# Nominal alginate:chitosan:calcium composition of each nanoparticle model.
SYSTEMS = ["1-1", "1-2", "2-1"]
SYSTEM_LABEL = {"1-1": "1:1", "1-2": "1:2", "2-1": "2:1"}
CHAIN_COUNTS = {
    "1-1": {"alginate": 30, "chitosan": 30, "calcium": 120},
    "1-2": {"alginate": 30, "chitosan": 60, "calcium": 120},
    "2-1": {"alginate": 60, "chitosan": 30, "calcium": 240},
}
# Antibiotics simulated in each system (8 copies of each).
PANEL = {
    "1-1": ["ATC", "MTC", "CPZ"],
    "1-2": ["ATC", "MTC", "CPZ"],
    "2-1": ["ATC", "VGM", "VCM", "RFP"],
}
ANTIBIOTIC_NAME = {
    "ATC": "Anhydrotetracycline",
    "MTC": "Methacycline",
    "CPZ": "Cefoperazone",
    "VGM": "Virginiamycin M1",
    "VCM": "Vancomycin",
    "RFP": "Rifampicin",
}
REPLICATES = [1, 2, 3]
COPIES_PER_ANTIBIOTIC = 8

# --- Analysis parameters -----------------------------------------------------
CONTACT_CUTOFF_A = 3.5      # heavy-atom distance defining direct contact
FINAL_WINDOW_NS = 20.0      # trailing window used for equilibrated averages
TRAJECTORY_NS = 100.0
STABLE_BINDING_FRACTION = 0.5   # contact fraction above which a molecule is "stably bound"
LONG_EVENT_NS = 10.0        # threshold for a long-lived contact event

# Radius-referenced classifier settings, reproduced from the trajectory script.
# Retained only to regenerate the comparison in Table S5; not used to measure adsorption.
POLYMER_RADIUS_PERCENTILE = 0.95
CORE_FRACTION = 45.0 / 55.0

# --- Plotting ----------------------------------------------------------------
SYSTEM_COLOR = {"1-1": "#2c6fbb", "1-2": "#c8622a", "2-1": "#3a7d44"}
ANTIBIOTIC_COLOR = {
    "ATC": "#1f77b4", "MTC": "#8c564b", "CPZ": "#d62728",
    "VGM": "#9467bd", "VCM": "#2ca02c", "RFP": "#ff7f0e",
}
FIGURE_DPI = 300
