#!/usr/bin/env python3
"""List actual cell coverage in pinned source files; no downloads or model fitting.

--source-root points to a workspace containing the full data/raw and data/unpacked
directories. The small output tables are written to this repository's data/catalog.
Names and identifiers are preserved; no spelling-based identity merges are made.
"""
from pathlib import Path
import argparse
import csv
from collections import Counter
import gzip
import json
import zipfile
import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / 'data/catalog'
SUMMARY = {}


def save(name, table, **facts):
    assert not table.empty, name
    table = table.sort_values('cell_name', kind='stable').reset_index(drop=True)
    table.to_csv(OUTPUT / f'{name}.csv', index=False)
    SUMMARY[name] = {'listed_rows':len(table), **facts}
    print(name, json.dumps(SUMMARY[name], ensure_ascii=False), flush=True)


def column(group, key):
    item = group[key]
    if isinstance(item, h5py.Dataset):
        return item.asstr()[:] if item.dtype.kind in 'OSU' else item[:]
    cats = item['categories']
    cats = cats.asstr()[:] if cats.dtype.kind in 'OSU' else cats[:]
    codes = item['codes'][:]
    result = np.full(len(codes), None, dtype=object)
    valid = codes >= 0
    result[valid] = cats[codes[valid]]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    args = parser.parse_args()
    source = args.source_root.resolve() / 'data'
    OUTPUT.mkdir(parents=True, exist_ok=True)
    model = pd.read_csv(source/'raw/CCLE_DepMap_24Q4/Model.csv', dtype=str).set_index('ModelID')

    def depmap(ids):
        return model.reindex(ids)[['CellLineName','OncotreeLineage','RRID']].rename(
            columns={'CellLineName':'cell_name','OncotreeLineage':'tissue','RRID':'cellosaurus_id'}).rename_axis('model_id').reset_index()

    expression = source/'raw/CCLE_DepMap_24Q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv'
    ids = pd.read_csv(expression,usecols=[0]).iloc[:,0].astype(str).tolist()
    table = depmap(ids)
    assert table.cell_name.notna().all() and len(set(ids)) == len(ids)
    table['expression_profiles'] = 1
    save('ccle_24q4',table,models_with_RNA=len(ids),models_in_Model_table=len(model),source_file=str(expression.relative_to(source)))

    for phase in ['LINCS_GSE70138','LINCS_GSE92742']:
        folder = source/'raw'/phase
        annotation = pd.read_csv(next(folder.glob('*cell_info*')),sep='\t',dtype=str).set_index('cell_id')
        signatures = pd.read_csv(next(folder.glob('*sig_info*')),sep='\t',usecols=['cell_id','pert_type'])
        counts = signatures.groupby(['cell_id','pert_type']).size().unstack(fill_value=0)
        table = pd.DataFrame({'cell_id':counts.index,'cell_name':counts.index,
                              'all_signatures':counts.sum(axis=1).to_numpy(),
                              'small_molecule_signatures':counts.get('trt_cp',pd.Series(0,index=counts.index)).to_numpy()})
        table['primary_site'] = table.cell_id.map(annotation.primary_site)
        table['source_cell_type'] = table.cell_id.map(annotation.cell_type)
        table['sample_type'] = table.cell_id.map(annotation.sample_type)
        save(phase.lower(),table,all_cells=signatures.cell_id.nunique(),
             cells_with_small_molecule_signatures=int(table.small_molecule_signatures.gt(0).sum()),
             all_signatures=len(signatures),small_molecule_signatures=int(table.small_molecule_signatures.sum()),
             source_file=str(next(folder.glob('*sig_info*')).relative_to(source)))

    bio = pd.read_csv(ROOT/'data/reference/hms_biorepeats.csv')
    conditions = pd.read_csv(ROOT/'data/reference/hms_conditions.csv')
    table = bio.groupby('cell_id').agg(raw_well_rows=('n_technical_rows','sum'),biorepeat_records=('source_file','size')).reset_index()
    table['cell_name'] = table.cell_id
    table['library_conditions'] = table.cell_id.map(conditions.groupby('cell_id').library_condition_id.nunique())
    table['distinct_cell_drug_dose'] = table.cell_id.map(conditions.groupby('cell_id').biological_condition_group.nunique())
    table['model_id'] = table.cell_id.map(conditions.drop_duplicates('cell_id').set_index('cell_id').ModelID)
    training = pd.read_csv(ROOT/'data/training/samples.csv')
    table['current_training_observations'] = table.cell_id.map(training.groupby('cell_id').size()).fillna(0).astype(int)
    save('hms',table,raw_wells=int(table.raw_well_rows.sum()),cells=len(table),source_file='HMS three Excel files; verified library/biorepeat references')

    folder = source/'unpacked/CTRPv2_2015'
    meta = pd.read_csv(folder/'v20.meta.per_cell_line.txt',sep='\t',dtype=str)
    experiments = pd.read_csv(folder/'v20.meta.per_experiment.txt',sep='\t')
    assert experiments.groupby('experiment_id').master_ccl_id.nunique().le(1).all()
    counts = Counter()
    for chunk in pd.read_csv(folder/'v20.data.per_cpd_post_qc.txt',sep='\t',usecols=['experiment_id'],chunksize=500000):
        counts.update({int(k):int(v) for k,v in chunk.experiment_id.value_counts().items()})
    pairs = pd.DataFrame({'experiment_id':list(counts),'dose_rows':list(counts.values())}).merge(
        experiments[['experiment_id','master_ccl_id']].drop_duplicates(),on='experiment_id',validate='one_to_one')
    table = pairs.groupby('master_ccl_id').dose_rows.sum().reset_index()
    table.master_ccl_id = table.master_ccl_id.astype(str)
    table = table.merge(meta[['master_ccl_id','ccl_name','ccle_primary_site']],on='master_ccl_id',validate='one_to_one')
    table = table.rename(columns={'ccl_name':'cell_name','ccle_primary_site':'tissue'})
    save('ctrpv2_2015',table,observed_cells=len(table),cell_metadata_rows=len(meta),dose_rows=int(table.dose_rows.sum()),source_file='v20.data.per_cpd_post_qc.txt + experiment/cell metadata')

    folder = source/'unpacked/CTRPv2_2016'
    meta = pd.read_csv(folder/'v21.meta.per_cell_line.txt',sep='\t',dtype=str)
    covered = {}
    for name, filename in [('auc_rows','v21.data.auc_sensitivities.txt'),
                           ('baseline_expression_features','v21.data.gex_avg_log2.txt'),
                           ('copy_number_features','v21.data.cnv_avg_log2.txt')]:
        counts = Counter()
        for chunk in pd.read_csv(folder/filename,sep='\t',usecols=['master_ccl_id'],dtype=str,chunksize=500000):
            counts.update({k:int(v) for k,v in chunk.master_ccl_id.value_counts().items()})
        covered[name] = pd.Series(counts,name=name)
    table = pd.concat(covered.values(),axis=1).fillna(0).astype(int).rename_axis('master_ccl_id').reset_index()
    table = table.merge(meta[['master_ccl_id','ccl_name','ccle_primary_site']],on='master_ccl_id',how='left',validate='one_to_one')
    table = table.rename(columns={'ccl_name':'cell_name','ccle_primary_site':'tissue'})
    assert table.cell_name.notna().all()
    auc = pd.read_csv(folder/'v21.data.auc_sensitivities.txt',sep='\t',usecols=['area_under_curve'])
    save('ctrpv2_2016',table,auc_cells=int(table.auc_rows.gt(0).sum()),
         expression_cells=int(table.baseline_expression_features.gt(0).sum()),
         copy_number_cells=int(table.copy_number_features.gt(0).sum()),
         auc_rows=int(table.auc_rows.sum()),auc_min=float(auc.area_under_curve.min()),auc_max=float(auc.area_under_curve.max()),
         source_file='CTRPv2 2016 AUC, baseline GEX and CNV tables; master_ccl_id metadata')

    for version in ['GDSC1','GDSC2']:
        path = next((source/'unpacked'/version).rglob('*.csv'))
        counts = Counter()
        for chunk in pd.read_csv(path,usecols=['CELL_LINE_NAME','SANGER_MODEL_ID'],dtype=str,chunksize=500000):
            counts.update({k:int(v) for k,v in chunk.fillna('').groupby(['CELL_LINE_NAME','SANGER_MODEL_ID']).size().items()})
        table = pd.DataFrame([{'cell_name':k[0],'sanger_model_id':k[1],'raw_well_rows':v} for k,v in counts.items()])
        save(version.lower(),table,distinct_raw_names=table.cell_name.nunique(),distinct_nonempty_Sanger_IDs=table.loc[table.sanger_model_id.ne(''),'sanger_model_id'].nunique(),raw_wells=int(table.raw_well_rows.sum()),source_file=str(path.relative_to(source)))

    for screen in ['primary','secondary']:
        folder = source/'raw/PRISM_19Q4'
        matrix = pd.read_csv(folder/f'{screen}-screen-replicate-collapsed-logfold-change.csv',index_col=0)
        info = pd.read_csv(folder/f'{screen}-screen-cell-line-info.csv',dtype=str).set_index('row_name')
        table = info.reindex(matrix.index)[['depmap_id','ccle_name','primary_tissue','passed_str_profiling']].reset_index()
        table.columns = ['matrix_row_id','model_id','cell_name','tissue','passed_str_profiling']
        table['nonmissing_measurements'] = matrix.notna().sum(axis=1).to_numpy()
        table['screen_metadata_row_present'] = table.matrix_row_id.isin(info.index)
        primary_info = pd.read_csv(folder/'primary-screen-cell-line-info.csv',dtype=str).set_index('row_name')
        table['originally_assigned_name'] = table.matrix_row_id.map(primary_info.originally_assigned_ccle_name)
        table['STR_notes_from_primary'] = table.matrix_row_id.map(primary_info.str_profiling_notes)
        table['failed_STR_label'] = table.matrix_row_id.str.endswith('_FAILED_STR')
        assert (table.cell_name.notna() | table.failed_STR_label).all()
        # Preserve failed-identity labels as unresolved assay rows. An originally
        # assigned name is not promoted to a verified current cell identity.
        named = int(table.cell_name.notna().sum())
        table.cell_name = table.cell_name.fillna(table.matrix_row_id)
        save(f'prism_{screen}',table,actual_matrix_rows=len(matrix),matrix_columns=len(matrix.columns),
             named_current_cell_rows=named,failed_STR_rows=int(table.failed_STR_label.sum()),
             valid_measurements=int(table.nonmissing_measurements.sum()),cell_metadata_rows=len(info),source_file=f'{screen}-screen-replicate-collapsed-logfold-change.csv')

    counts = Counter()
    barcode_count = 0
    for archive in sorted((source/'raw/MIXseq').glob('*.zip')):
        with zipfile.ZipFile(archive) as z:
            for filename in z.namelist():
                if not filename.endswith('/classifications.csv') or filename.startswith('__MACOSX/'):
                    continue
                prefix = filename.rsplit('/',1)[0]+'/'
                metadata = pd.read_csv(z.open(filename))
                barcodes = pd.read_csv(z.open(prefix+'barcodes.tsv'),header=None,sep='\t')[0]
                aligned = metadata.set_index('barcode').reindex(barcodes)
                assert aligned.cell_quality.notna().all()
                barcode_count += len(barcodes)
                selected = aligned[aligned.cell_quality.eq('normal') & aligned.DepMap_ID.notna()]
                counts.update({k:int(v) for k,v in selected.DepMap_ID.value_counts().items()})
    counts = pd.Series(counts).sort_index()
    table = depmap(counts.index.tolist())
    table['normal_cells_with_ID'] = counts.to_numpy()
    # An ID without a name in this release remains visible under its source ID.
    table.cell_name = table.cell_name.fillna(table.model_id)
    table['identity_name_source'] = 'DepMap24Q4 Model.csv matched by exact DepMap ID'
    save('mixseq',table,matrix_barcodes=barcode_count,normal_cells=int(counts.sum()),actual_DepMap_IDs=len(counts),source_file='26 author archives/classifications.csv; author normal flag and nonmissing DepMap_ID')

    for number in [2,3]:
        path = source/f'raw/sciPlex/SrivatsanTrapnell2020_sciplex{number}.h5ad'
        with h5py.File(path) as f:
            names = pd.Series(column(f['obs'],'cell_line'))
            counts = names.value_counts()
            table = pd.DataFrame({'cell_name':counts.index,'cell_or_nucleus_rows':counts.to_numpy()})
            facts = {'annotated_rows':int(counts.sum()),'missing_cell_line_rows':int(names.isna().sum()),'total_rows':len(names)}
            if number==2:
                doses = pd.Series(column(f['obs'],'dose_value'))
                tags = pd.Series(column(f['obs'],'top_oligo'))
                derived = tags.str.split('_').str[1]
                assert (derived[doses.notna()] == doses[doses.notna()]).all()
                facts.update(raw_dose_labels=sorted(doses.dropna().unique().tolist()),missing_dose_labels=int(doses.isna().sum()),dose_unit_field_present='dose_unit' in f['obs'],nonmissing_dose_labels_match_top_oligo_split=True)
            if number==3:
                times = pd.Series(column(f['obs'],'time'))
                facts['cell_time_counts'] = pd.DataFrame({'cell':names,'time':times}).groupby(['cell','time'],dropna=False).size().rename('rows').reset_index().astype(object).where(lambda x:pd.notna(x),None).to_dict('records')
        save(f'sciplex{number}',table,**facts,source_file=path.name)

    header = {}
    with gzip.open(source/'raw/NCI_TPW/GSE116436_series_matrix.txt.gz','rt') as stream:
        for line in stream:
            fields = next(csv.reader([line],delimiter='\t'))
            if fields and fields[0] in ['!Sample_title','!Sample_geo_accession','!Sample_source_name_ch1']:
                header[fields[0].removeprefix('!Sample_')] = fields[1:]
            if fields and fields[0] == '!series_matrix_table_begin':
                break
    nci = pd.DataFrame(header)
    nci['dose_nM'] = pd.to_numeric(nci.title.str.extract(r'_([\d.]+)nM_[\d.]+h$')[0],errors='coerce')
    assert nci.dose_nM.notna().all()
    table = nci.groupby('source_name_ch1').agg(expression_samples=('geo_accession','size'),zero_dose_samples=('dose_nM',lambda x:int(x.eq(0).sum()))).reset_index()
    table['cell_name'] = table.source_name_ch1.str.removeprefix('human tumor cell line ')
    save('nci_tpw',table,expression_samples=len(nci),cells=len(table),source_file='GSE116436_series_matrix.txt.gz sample headers')
    save('drugseq',pd.DataFrame([{'cell_name':'U2OS','profiled_well_columns':59904,'processed_sample_members':56628,'cell_identity_source':'Novartis author README; not inferred from sample count'}]),cells=1,source_file='Exp_gzip.RData + comparisons_metadata.RData; cell identity from Novartis README')

    folder = source/'raw/Tahoe_100M'
    annotation = pq.read_table(folder/'cell_line_metadata.parquet').to_pandas()
    identity = annotation[['cell_name','Cell_ID_DepMap','Cell_ID_Cellosaur','Organ']].drop_duplicates()
    assert identity.Cell_ID_Cellosaur.is_unique
    shard = pq.read_table(folder/'train-00000-of-03388.parquet',columns=['cell_line_id','sample','drug']).to_pandas()
    counts = shard.cell_line_id.value_counts().rename_axis('Cell_ID_Cellosaur').rename('observed_shard_cells').reset_index()
    table = counts.merge(identity,on='Cell_ID_Cellosaur',how='left',validate='one_to_one')
    assert table.cell_name.notna().all()
    save('tahoe_shard',table,metadata_rows=len(annotation),metadata_unique_identities=len(identity),observed_cell_identities=len(table),shard_rows=len(shard),observed_samples=shard['sample'].nunique(),observed_treatment_names=shard.drug.nunique(),source_file='train-00000-of-03388.parquet + cell_line_metadata.parquet')
    put = {'date':'2026-11-09','counts_are':'actual cells in the specified files; no alias collapse','datasets':SUMMARY}
    (OUTPUT/'summary.json').write_text(json.dumps(put,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
