#!/usr/bin/env python3
"""
Convert source spatial transcriptomics data into app-ready CSV files.

Usage (set VGEP_DATA to the data/ folder of the Zenodo archive, 10.5281/zenodo.19633285):
    python scripts/convert_data.py

Outputs (in data/):
    Atlas.csv          — 24 training genes + coordX, coordY
    Validation.csv     — 23 validation genes (7 protein + 13 smiFISH + Tollo, yrt, stan) + coordX, coordY
    vGEP.csv           — 4,312 genes + coordX, coordY
    macrochaete.csv    — 8 macrochaetae positions (coordX;coordY)
    md_domains.csv     — 5-MD morphodynamic domain (1-5) per grid bin + coordX, coordY
"""

import os
import pandas as pd
import mudata

# ============================================================
# Configuration — edit paths here to point to new data sources
# ============================================================

DATA_ROOT = os.environ.get(
    "VGEP_DATA", os.path.join(os.path.dirname(__file__), "..", "..", "..", "data")
)
ATLAS_VALIDATION_PATH = os.path.join(DATA_ROOT, "morphology", "atlas_validation.txt")
H5MU_PATH = os.path.join(
    DATA_ROOT, "gep_results", "spatial_transcriptome_output", "best_transcriptome.h5mu"
)
MACROCHAETE_PATH = os.path.join(DATA_ROOT, "morphology", "macrochaete.csv")
MD_PATH = os.path.join(
    DATA_ROOT, "..", "figures", "fig_morphodynamic_domains", "clustering", "cluster_labels_5.csv"
)

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "data")

# Column renames: file column name → canonical gene name
RENAME = {
    "BH1": "B-H1",
    "Dpp": "dpp",
    "robo": "robo1",
    "slit": "sli",
    "gene18w": "18w",
}

TRAINING_GENES = [
    "al", "ara", "B-H1", "bi", "Diap1", "esg", "exd", "eyg",
    "h", "hh", "lgs", "mirr", "ptc", "salm", "sd", "tsh",
    "ush", "sr", "usp", "wg", "en", "rpr", "hid", "Pax",
]

VALIDATION_GENES = [
    "C15", "caup", "crp", "Dl", "dpp", "fj", "hbs", "hth",
    "nw", "pdm3", "Pvr", "Qsox2", "robo1", "sli", "slow",
    "toe", "Wnt6", "svp", "vvl", "18w", "Tollo", "yrt", "stan",
]


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ----------------------------------------------------------
    # 1. Atlas & Validation from atlas_validation.txt
    # ----------------------------------------------------------
    print("Reading atlas_validation.txt ...")
    av = pd.read_csv(ATLAS_VALIDATION_PATH, sep=",")
    av.rename(columns=RENAME, inplace=True)
    print(f"  Raw shape: {av.shape}")

    # Verify all expected genes are present
    for gene in TRAINING_GENES:
        assert gene in av.columns, f"Training gene '{gene}' not found in columns"
    for gene in VALIDATION_GENES:
        assert gene in av.columns, f"Validation gene '{gene}' not found in columns"

    # Mask: rows where ALL training genes are non-NaN
    mask = av[TRAINING_GENES].notna().all(axis=1)
    print(f"  Rows with all training genes non-NaN: {mask.sum()}")

    # Atlas.csv
    atlas = av.loc[mask, TRAINING_GENES + ["coordX", "coordY"]].reset_index(drop=True)
    atlas_path = os.path.join(OUTPUT_DIR, "Atlas.csv")
    atlas.to_csv(atlas_path, index=False)
    print(f"  Atlas.csv: {atlas.shape} -> {atlas_path}")

    # Validation.csv (same row mask; some validation genes may have NaN — that's OK)
    validation = av.loc[mask, VALIDATION_GENES + ["coordX", "coordY"]].reset_index(drop=True)
    val_path = os.path.join(OUTPUT_DIR, "Validation.csv")
    validation.to_csv(val_path, index=False)
    print(f"  Validation.csv: {validation.shape} -> {val_path}")

    # ----------------------------------------------------------
    # 2. vGEP from h5mu
    # ----------------------------------------------------------
    print("\nReading best_transcriptome.h5mu ...")
    mdata = mudata.read_h5mu(H5MU_PATH)
    st = mdata.mod["st"]
    print(f"  ST modality: {st.shape[0]} locations x {st.shape[1]} genes")

    # Expression matrix
    X = st.X
    if hasattr(X, "toarray"):
        X = X.toarray()

    # Gene names (strip "st:" prefix if present)
    gene_names = [g.replace("st:", "") for g in st.var_names]

    # Coordinates
    coords_x = st.obs["x"].astype(float).values
    coords_y = st.obs["y"].astype(float).values

    # Build DataFrame
    vgep = pd.DataFrame(X, columns=gene_names)
    vgep["coordX"] = coords_x
    vgep["coordY"] = coords_y

    vgep_path = os.path.join(OUTPUT_DIR, "vGEP.csv")
    vgep.to_csv(vgep_path, index=False)
    print(f"  vGEP.csv: {vgep.shape} -> {vgep_path}")

    # ----------------------------------------------------------
    # 3. Macrochaetae
    # ----------------------------------------------------------
    print("\nReading macrochaete.csv ...")
    macro = pd.read_csv(MACROCHAETE_PATH, header=None)
    print(f"  Raw shape: {macro.shape}")

    # Transpose: 2 rows x N cols → N rows x 2 cols
    # Row 0 = Y coords, Row 1 = X coords
    macro_t = pd.DataFrame({
        "coordX": macro.iloc[1].values,
        "coordY": macro.iloc[0].values,
    })
    macro_path = os.path.join(OUTPUT_DIR, "macrochaete.csv")
    macro_t.to_csv(macro_path, index=False, sep=";")
    print(f"  macrochaete.csv: {macro_t.shape} -> {macro_path}")

    # ----------------------------------------------------------
    # 4. Morphodynamic domains (cluster id 0-4 = 5-MD1..5), restricted to the app grid
    # ----------------------------------------------------------
    print("\nReading cluster_labels_5.csv ...")
    md = pd.read_csv(MD_PATH)
    md.columns = ["coordX", "coordY", "domain"]
    md["domain"] += 1
    md = md.merge(vgep[["coordX", "coordY"]].astype(int), on=["coordX", "coordY"])
    md_path = os.path.join(OUTPUT_DIR, "md_domains.csv")
    md[["domain", "coordX", "coordY"]].to_csv(md_path, index=False)
    print(f"  md_domains.csv: {md.shape} -> {md_path}")

    # ----------------------------------------------------------
    # Summary
    # ----------------------------------------------------------
    print("\n=== Summary ===")
    print(f"  Atlas:       {atlas.shape[0]} rows x {atlas.shape[1]} cols "
          f"({atlas.shape[1] - 2} genes)")
    print(f"  Validation:  {validation.shape[0]} rows x {validation.shape[1]} cols "
          f"({validation.shape[1] - 2} genes)")
    print(f"  vGEP:        {vgep.shape[0]} rows x {vgep.shape[1]} cols "
          f"({vgep.shape[1] - 2} genes)")
    print(f"  Macrochaete: {macro_t.shape[0]} rows x {macro_t.shape[1]} cols")


if __name__ == "__main__":
    main()
