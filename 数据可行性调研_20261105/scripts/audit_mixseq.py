#!/usr/bin/env python3
"""Check four original MIX-seq archives and form population-level expression pairs."""
from pathlib import Path
import json
import zipfile
import xml.etree.ElementTree as ET
import pandas as pd
from scipy import sparse
from scipy.io import mmread

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "audit"


def main():
    # Read the original supplementary workbook without requiring an Excel package.
    book = zipfile.ZipFile(ROOT / "metadata/mixseq_Supplementary_Tables.xlsx")
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    shared = ["".join(x.itertext()) for x in ET.fromstring(book.read("xl/sharedStrings.xml")).findall("m:si", ns)]
    drug_rows = []
    for row in ET.fromstring(book.read("xl/worksheets/sheet3.xml")).findall(".//m:row", ns):
        values = []
        for cell in row.findall("m:c", ns):
            v = cell.find("m:v", ns)
            if v is not None:
                values.append(shared[int(v.text)] if cell.attrib.get("t") == "s" else v.text)
        if values and values[0] == "Trametinib":
            drug_rows.append(values)
    assert len(drug_rows) == 1, drug_rows
    trametinib_dose_uM = float(drug_rows[0][1])
    summaries, counts = [], []
    for drug in ["DMSO", "Trametinib"]:
        for time in [6, 24]:
            name = f"{drug}_{time}hr_expt1"
            z = zipfile.ZipFile(ROOT / "metadata" / f"mixseq_{name}.zip")
            prefix = name + "/"
            meta = pd.read_csv(z.open(prefix + "classifications.csv"))
            barcodes = pd.read_csv(z.open(prefix + "barcodes.tsv"), sep="\t", header=None)[0]
            genes = pd.read_csv(z.open(prefix + "genes.tsv"), sep="\t", header=None)
            m = mmread(z.open(prefix + "matrix.mtx")).tocsr()
            assert m.shape == (len(genes), len(barcodes)), (name, m.shape)
            assert barcodes.is_unique and meta.barcode.is_unique, name
            aligned = meta.set_index("barcode").reindex(barcodes)
            assert aligned.cell_quality.notna().all(), name
            valid = aligned.cell_quality.eq("normal") & aligned.DepMap_ID.notna()
            good = aligned[valid]
            n = good.groupby("DepMap_ID").size()
            n = n.rename("n_cells").reset_index()
            n["condition"] = name
            n["drug"] = drug
            n["time_h"] = time
            n["expression_dose_uM"] = trametinib_dose_uM if drug == "Trametinib" else 0.
            counts.append(n)
            summaries.append({"condition": name, "matrix_genes": m.shape[0],
                              "matrix_cells": m.shape[1], "qc_normal_cells": int(valid.sum()),
                              "qc_cell_lines": len(n)})
            # Sum raw counts per cell line. Cell-level expression is never a viability label.
            ids = sorted(good.DepMap_ID.unique())
            lookup = {x: i for i, x in enumerate(ids)}
            positions = [i for i, ok in enumerate(valid) if ok]
            grouping = sparse.csr_matrix(([1] * len(positions),
                                          ([lookup[aligned.iloc[i].DepMap_ID] for i in positions], positions)),
                                         shape=(len(ids), len(barcodes)))
            pseudo = grouping @ m.T
            sparse.save_npz(OUT / f"mixseq_{name}_pseudobulk_counts.npz", pseudo)
            pd.Series(ids, name="DepMap_ID").to_csv(OUT / f"mixseq_{name}_pseudobulk_cells.csv", index=False)
            genes.to_csv(OUT / f"mixseq_{name}_genes.tsv", sep="\t", index=False, header=False)
    all_counts = pd.concat(counts, ignore_index=True)
    all_counts.to_csv(OUT / "mixseq_pilot_cell_counts.csv", index=False)
    wide = all_counts.pivot(index="DepMap_ID", columns="condition", values="n_cells").fillna(0).astype(int)
    lab = pd.read_csv(OUT / "mixseq_feature_labels.csv")
    lab = lab[lab.comparison.eq("Trametinib_24hr_expt1") & lab.in_pool.eq(True)]
    assert not lab.DEPMAP_ID.duplicated().any()
    pairs = wide.merge(lab[["DEPMAP_ID", "sens"]], left_index=True, right_on="DEPMAP_ID", how="left").rename(columns={"DEPMAP_ID": "DepMap_ID"})
    pairs["all_four_conditions_ge_10_cells"] = pairs[wide.columns].ge(10).all(axis=1)
    pairs["all_four_conditions_ge_20_cells"] = pairs[wide.columns].ge(20).all(axis=1)
    pairs["sensitivity_label_type"] = "1_minus_quantile_normalized_cross_screen_AUC; not_dose_specific_viability"
    pairs.to_csv(OUT / "mixseq_pilot_expression_sensitivity_pairs.csv", index=False)
    result = {"archives": summaries, "trametinib_expression_dose_uM": trametinib_dose_uM,
              "dose_source": "Original Figshare Supplementary Tables.xlsx, sheet3, Trametinib row",
              "qc_normal_total": sum(x["qc_normal_cells"] for x in summaries),
              "cell_lines_in_any_condition": len(wide),
              "all_four_conditions_ge_10_cells": int(pairs.all_four_conditions_ge_10_cells.sum()),
              "all_four_conditions_ge_20_cells": int(pairs.all_four_conditions_ge_20_cells.sum()),
              "ge_20_and_sensitivity_label": int((pairs.all_four_conditions_ge_20_cells & pairs.sens.notna()).sum()),
              "note": "One drug, two destructive expression sampling times, time-matched controls. Population-level conditions; no individual-cell trajectories or concentration-specific viability labels."}
    (OUT / "mixseq_pilot_audit.json").write_text(json.dumps(result, ensure_ascii=False, indent=2))
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
