#!/usr/bin/env python3
"""Audit downloaded source tables. Counts describe metadata, not expression QC."""
from pathlib import Path
import hashlib
import json
import re
import sys
import zipfile

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "metadata"
OUT = ROOT / "audit"
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / "evidence/xlrd-2.0.2-py2.py3-none-any.whl"))


def dump(name, obj):
    (OUT / name).write_text(json.dumps(obj, ensure_ascii=False, indent=2,
                                     default=lambda x: x.item() if isinstance(x, np.generic) else str(x)), encoding="utf-8")


def schema(name, df):
    return [{"source_file": name, "field": str(c), "dtype": str(df[c].dtype),
             "rows": len(df), "non_null": int(df[c].notna().sum())} for c in df.columns]


def lincs(path, source):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    fields = schema(path.name, df)
    df = df[df.pert_type.eq("trt_cp")].copy()
    # Phase I has quantitative doses; pert_idose is a rounded display label.
    if "pert_dose" in df:
        df["dose_uM"] = pd.to_numeric(df.pert_dose, errors="coerce")
        units = df.pert_dose_unit.astype(str).str.replace("µ", "u").str.replace("μ", "u").str.lower()
        factors = units.map({"um": 1., "nm": .001, "mm": 1000., "m": 1e6})
        df["dose_uM"] *= factors
    else:
        dose = df.pert_idose.astype(str).str.extract(r"^\s*([\d.eE+-]+)\s*([^\s]+)")
        units = dose[1].str.replace("µ", "u").str.replace("μ", "u").str.lower()
        df["dose_uM"] = pd.to_numeric(dose[0], errors="coerce") * units.map({"um": 1., "nm": .001, "mm": 1000., "m": 1e6})
    df["time_h"] = pd.to_numeric(df.pert_itime.str.extract(r"^([\d.]+)")[0], errors="coerce")
    df["dose_uM"] = df.dose_uM.round(6)
    df["source"] = source
    return df[df.dose_uM.gt(0) & df.time_h.gt(0)].copy(), fields


def main():
    result = {"audit_date": "2026-11-05", "scope": "Downloaded metadata and supplementary tables; LINCS expression matrices not downloaded or quality-filtered. See separate MIX-seq audit for original expression matrix checks."}
    fields = []
    frames = []
    for source in ["GSE92742", "GSE70138"]:
        df, f = lincs(META / f"{source}_sig_info.txt.gz", source)
        fields += f
        frames.append(df)
        result[source] = {"positive_dose_chemical_signatures": len(df),
                          "cells": df.cell_id.nunique(), "drugs": df.pert_id.nunique(),
                          "time_counts": {str(k): int(v) for k, v in df.time_h.value_counts().sort_index().items()}}
    lin = pd.concat(frames, ignore_index=True)
    ljp = lin[lin.sig_id.str.startswith("LJP")].copy()
    hfile = "Screen20237_LJP_processedData1.xls"
    h = pd.read_excel(META / hfile, sheet_name="Data")
    fields += schema(hfile + ":Data", h)
    h.to_csv(OUT / "hms_screen20237_observations.csv", index=False)
    h["dose_uM"] = pd.to_numeric(h["Small Mol Concentration"]).round(6)
    h = h.rename(columns={"CellLine": "cell_id", "BRDid": "pert_id"})
    keys = ["cell_id", "pert_id", "dose_uM"]
    grouped = h.groupby(keys, dropna=False).agg(
        n_observations=("Relative Cell Count", "size"),
        relative_cell_count=("Relative Cell Count", "mean"),
        growth_rate_value=("Relative Growth Rate After Treatment", "mean"),
        viability_time_h=("Timepoint", "first"),
        drug_name=("DrugName", "first"),
    ).reset_index()
    grouped.to_csv(OUT / "hms_screen20237_conditions.csv", index=False)
    matched = ljp.merge(grouped, on=keys, how="inner", validate="many_to_one")
    matched["pairing_status"] = "coordinated_project_same_cell_drug_dose_parallel_assays_provenance_to_check"
    matched.to_csv(OUT / "hms_lincs_candidate_pairs.csv", index=False)
    unique = matched.drop_duplicates(keys + ["time_h"])
    # A cross-time candidate must be within one release, not assembled across releases.
    time_sets = matched.groupby(keys + ["source"]).time_h.agg(lambda x: set(x))
    times = sorted(unique.time_h.unique())
    result["HMS_LINCS"] = {
        "raw_observations": len(h), "cells": h.cell_id.nunique(),
        "cell_names": sorted(h.cell_id.unique()), "compounds": h.pert_id.nunique(),
        "dose_values_uM": sorted(h.dose_uM.unique()), "functional_endpoint_h": sorted(h.Timepoint.unique()),
        "unique_functional_conditions": len(grouped),
        "matched_LJP_signature_records": len(matched),
        "matched_unique_expression_conditions": len(unique),
        "matched_unique_functional_conditions": len(unique.drop_duplicates(keys)),
        "matched_cell_drug_pairs": len(unique.drop_duplicates(keys[:2])),
        "expression_times_h": times,
        "by_source_and_time": matched.groupby(["source", "time_h"]).size().rename("n").reset_index().to_dict("records"),
        "unique_conditions_by_source_and_time": matched.drop_duplicates(keys + ["source", "time_h"]).groupby(["source", "time_h"]).size().rename("n").reset_index().to_dict("records"),
        "matched_3_to_24_conditions": int(time_sets.apply(lambda s: {3., 24.}.issubset(s)).sum()),
        "matched_6_to_24_conditions": int(time_sets.apply(lambda s: {6., 24.}.issubset(s)).sum()),
        "qualification": "Condition correspondence from IDs, doses and LJP prefix. Different wells/timepoints; no gene-matrix QC or batch-level provenance confirmation yet. Do not sum rows as independent trajectories.",
    }
    drug = h[["pert_id", "drug_name"]].copy() if "drug_name" in h else h[["pert_id", "DrugName"]].rename(columns={"DrugName": "drug_name"})
    drug = drug.drop_duplicates()
    pert = pd.read_csv(META / "GSE92742_pert_info.txt.gz", sep="\t", low_memory=False)
    mapping = drug.merge(pert.drop_duplicates("pert_id"), on="pert_id", how="left", validate="one_to_one")
    mapping.to_csv(OUT / "hms_drug_mapping.csv", index=False)
    models = pd.read_csv(META / "depmap24q4_Model.csv")
    defaults = pd.read_csv(META / "depmap24q4_OmicsDefaultModelProfiles.csv")
    candidate_models = models[models.StrippedCellLineName.isin(h.cell_id.unique())]
    cells = pd.DataFrame({"cell_id": sorted(h.cell_id.unique())}).merge(
        candidate_models[["ModelID", "CellLineName", "StrippedCellLineName", "RRID"]],
        left_on="cell_id", right_on="StrippedCellLineName", how="left", validate="one_to_one")
    rna = defaults[defaults.ProfileType.eq("rna")]
    cells = cells.merge(rna[["ModelID", "ProfileID"]], on="ModelID", how="left", validate="one_to_one")
    cells.to_csv(OUT / "hms_cell_mapping.csv", index=False)
    structures_present = mapping.canonical_smiles.notna() & ~mapping.canonical_smiles.astype(str).isin(["-666", "", "nan"])
    complete = matched[matched.pert_id.isin(mapping.loc[structures_present, "pert_id"]) &
                       matched.cell_id.isin(cells.loc[cells.ProfileID.notna(), "cell_id"])].copy()
    complete.to_csv(OUT / "hms_complete_input_metadata_candidate_pairs.csv", index=False)
    complete_times = complete.groupby(keys + ["source"]).time_h.agg(set)
    result["HMS_input_metadata_availability"] = {
        "compounds_with_structure_string_in_downloaded_metadata": int(structures_present.sum()),
        "compounds_without_structure_string": int((~structures_present).sum()),
        "cells_with_default_rna_profile_index": int(cells.ProfileID.notna().sum()),
        "functional_conditions_with_both_input_metadata": len(complete.drop_duplicates(keys)),
        "within_release_3_to_24_conditions_with_both_input_metadata": int(complete_times.apply(lambda s: {3., 24.}.issubset(s)).sum()),
        "qualification": "Metadata completeness only, using Phase I structure table and pinned DepMap 24Q4 default RNA index. Missing structures may be recoverable elsewhere. Neither SMILES validity nor RNA/expression matrices have been checked."}
    result["HMS_DepMap_baseline_metadata"] = {"mapped_models": int(cells.ModelID.notna().sum()),
                                             "default_rna_profiles": int(cells.ProfileID.notna().sum()),
                                             "qualification": "Profile index availability only; RNA expression matrix not downloaded. Baseline from a different experiment is not a time-matched L1000 control."}
    fields += schema("depmap24q4_Model.csv", models) + schema("depmap24q4_OmicsDefaultModelProfiles.csv", defaults)
    # Paper's response-class tables are a separate cohort from Screen20237.
    supplement = META / "niepel2017_supplementary.zip"
    z = zipfile.ZipFile(supplement)
    paper_summary = {}
    for num in [5, 6, 7]:
        name = f"41467_2017_1383_MOESM{num}_ESM.txt"
        df = pd.read_csv(z.open(name), sep="\t", encoding="latin1")
        fields += schema(name, df)
        df.to_csv(OUT / f"niepel_supp{num}.csv", index=False)
        entry = {"rows": len(df), "columns": list(df.columns)}
        if "Response class" in df:
            entry.update(cells=df["Cell line"].nunique(), drugs=df["Drug Name"].nunique(),
                         unique_conditions=len(df.drop_duplicates(["Cell line", "Drug Name", "Concentration"])),
                         response_classes={str(k): int(v) for k, v in df["Response class"].value_counts(dropna=False).items()},
                         spearman_SCS_GR=float(df[["SCS", "GR value (72h)"]].corr(method="spearman").iloc[0, 1]))
            quadrants = pd.crosstab(df.SCS.gt(1.3).rename("SCS_gt_1_3"), df["GR value (72h)"].lt(.66).rename("GR_lt_0_66"))
            quadrants.to_csv(OUT / f"niepel_supp{num}_quadrants.csv")
            entry["quadrants"] = [{"SCS_gt_1_3": bool(i), "GR_lt_0_66": bool(j), "n": int(quadrants.loc[i, j])} for i in quadrants.index for j in quadrants.columns]
        paper_summary[f"supp{num}"] = entry
    result["Niepel_paper_supplements"] = paper_summary
    # PRISM metadata describe the tested design, not the measured nonmissing matrix.
    for stem in ["prism_secondary_treatment", "prism_secondary_cell"]:
        f = META / f"{stem}.csv"
        if f.exists():
            df = pd.read_csv(f, low_memory=False)
            fields += schema(f.name, df)
            result[stem] = {"rows": len(df), "columns": list(df.columns)}
            for c in ["name", "broad_id", "depmap_id", "dose", "ccle_name"]:
                if c in df:
                    result[stem][f"unique_{c}"] = df[c].nunique()
    for kind in ["cell_line", "drug", "sample"]:
        f = META / f"tahoe_{kind}_metadata.parquet"
        if f.exists() and f.read_bytes()[:4] == b"PAR1":
            df = pd.read_parquet(f)
            fields += schema(f.name, df)
            result[f"Tahoe_{kind}"] = {"rows": len(df), "columns": list(df.columns)}
            if kind == "cell_line":
                result[f"Tahoe_{kind}"].update(unique_cellosaurus=df.Cell_ID_Cellosaur.nunique(),
                                                unique_depmap=df.Cell_ID_DepMap.nunique(),
                                                qualification="Driver annotation rows include a broader panel. These IDs are not evidence of expression coverage; filter by IDs actually present in expression/obs data.")
            df.to_csv(OUT / f"tahoe_{kind}_metadata.csv", index=False)
    # All evidence files are inventoried with hashes, including access-gate HTML.
    manifest = []
    for folder in ["metadata", "evidence"]:
        for f in sorted((ROOT / folder).iterdir()):
            if not f.is_file():
                continue
            raw = f.read_bytes()
            manifest.append({"path": str(f.relative_to(ROOT)), "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(),
                             "html_response": raw.lstrip().lower().startswith((b"<!doctype html", b"<html"))})
    pd.DataFrame(manifest).to_csv(OUT / "downloaded_file_manifest.csv", index=False)
    pd.DataFrame(fields).to_csv(OUT / "field_dictionary.csv", index=False)
    dump("audit_results.json", result)
    print((OUT / "audit_results.json").read_text())


if __name__ == "__main__":
    main()
