#!/usr/bin/env python3
"""Build project-wide read-only-source inventories and decision registers."""
from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / '项目总设计_20261106'


def main():
    OUT.mkdir(exist_ok=True)
    status = json.loads((ROOT / 'data/download_status.json').read_text())
    records = []
    for r in status['files']:
        path = ROOT / 'data' / r['path']
        size = path.stat().st_size if path.is_file() else None
        recorded = r.get('bytes', r.get('downloaded_bytes'))
        records.append({'dataset': r['dataset'], 'name': r['name'], 'source_url': r['url'],
                        'path': str(path.relative_to(ROOT)), 'scope': r.get('scope'),
                        'recorded_status': r['status'], 'exists': path.is_file(),
                        'actual_bytes': size, 'recorded_bytes': recorded,
                        'size_matches_record': size == recorded if recorded is not None else None,
                        'checksum': r.get('checksum'), 'checksum_algorithm': r.get('checksum_algorithm'),
                        'official_checksum_verified_in_download_record': r.get('official_checksum_verified', False)})
    inventory = pd.DataFrame(records)
    assert len(inventory) == 78 and inventory.exists.all() and inventory.recorded_status.eq('已完成').all()
    assert inventory.size_matches_record.dropna().all()
    inventory.to_csv(OUT / 'file_inventory.csv', index=False)
    columns = ['dataset_id', 'source', 'source_url', 'local_path', 'observation_unit', 'readout_and_scale', 'conditions_and_controls', 'project_role', 'actual_stage', 'next_output']
    rows = [
        ['CCLE_24Q4', 'Broad DepMap 24Q4', 'https://depmap.org/portal/download/all/', 'data/raw/CCLE_DepMap_24Q4', '未处理细胞系背景', 'log2(TPM+1)', '独立研究背景；非HMS配对t0', '背景输入和未见细胞', '五细胞978基因提取完成', '身份/培养复核与扩展背景清单'],
        ['LINCS_I_II', 'Broad LINCS GEO', 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138', 'data/raw/LINCS_GSE70138|data/raw/LINCS_GSE92742', '实例/汇总签名', 'L3归一化表达；L5响应z', '细胞/药物/剂量/时间/板/vehicle', '主响应目标和额外表达监督', '两期下载完成；II已解压；HMS6814签名已提取；I大矩阵尚未数值提取', '按需要提取I；实例/vehicle QC和跨发布重复审计'],
        ['HMS', 'HMS LINCS Joint Project', 'https://lincs.hms.harvard.edu/db/datasets/20237/', 'data/raw/HMS_LINCS', '技术孔→库板内生物重复→条件', 'Nt/Nc；相对净增长', '72h；DMSO及N0', '主功能标签和表达功能关系', '原孔核查/汇总/v3数组及划分草案完成；生物QC未冻结', '入组与正式训练版本'],
        ['CTRPv2', 'Broad CTD2 2015扩展/2016补充', 'https://portals.broadinstitute.org/ctrp.v2/', 'data/raw/CTRPv2|data/unpacked/CTRPv2_2015', '实验曲线中的实测剂量点', 'cpd_avg_pv比例；拟合值/AUC分列', '72h ATP；DMSO', '扩展剂量功能监督', '6171005行扫描；最近剂量关系生成；全局去重/身份QC未完成', '保留原剂量的单药功能长表'],
        ['GDSC1', 'Sanger', 'https://www.cancerrxgene.org/', 'data/raw/GDSC', '板内原始孔及曲线', 'assay强度；拟合LN_IC50/AUC', '72h暴露；Resazurin/Syto60；controls', '来源特定功能监督', '原始5976183行和板级control检查完成', '官方规则归一化与小批实例核对'],
        ['GDSC2', 'Sanger', 'https://www.cancerrxgene.org/', 'data/raw/GDSC', '板内原始孔及曲线', 'CellTiter-Glo强度；拟合曲线', '72h暴露；流程Duration=4非96h给药', '来源特定功能监督', '原始15069449行检查完成', '官方规则归一化和身份/剂量对接'],
        ['PRISM_19Q4', 'Broad PRISM Repurposing', 'https://depmap.org/repurposing/', 'data/raw/PRISM_19Q4', '细胞系×药物剂量/筛选', '处理后的log2相对条形码丰度', '五天120h；DMSO；混合培养', '来源特定功能和迁移', 'secondary实际矩阵489x13008已检查', '身份STR/重做选择及长表导出'],
        ['MIXseq', 'McFarland等2020', 'https://doi.org/10.1038/s41467-020-17440-w', 'data/raw/MIXseq', 'barcode及实验重复伪bulk', 'scRNA counts；附带sens来自GDSC/PRISM', '实际实验药物/时间及DMSO；多时间点', '外部表达和时间/关联分解', '26个单药对照ZIP及3支持文件内容已核', '细胞QC/同批对照/伪bulk及锁定外部条件'],
        ['sciPlex2', 'Srivatsan等2020/scPerturb', 'https://doi.org/10.1126/science.aax6234', 'data/raw/sciPlex', '单细胞/核', 'scRNA count', 'A549四药24h及vehicle', '外部表达流程试验', '24262x58347 h5ad实际读取', '身份/基因/QC和重复级伪bulk'],
        ['sciPlex3', 'Srivatsan等2020/scPerturb', 'https://doi.org/10.1126/science.aax6234', 'data/raw/sciPlex', '单细胞/核', '人鼠特征counts；回收数量代理非全面功能真值', '三细胞；24/72h结构不平衡；nM剂量', '表达外部验证和异质性', '799317x110983元数据已核', '人基因/重复基因/缺注释/QC及条件锁定'],
        ['NCI_TPW', 'NCI GSE116436', 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116436', 'data/raw/NCI_TPW', '芯片样本', 'GPL571探针表达', '2/6/24h及零剂量对照', '外部表达与时间', '7623样本表头和条件已核', '探针映射/对照配对和稳定身份'],
        ['DRUGseq', 'Novartis U2OS MoABox', 'https://zenodo.org/records/14291446', 'data/raw/Novartis_DRUGseq', '药物处理孔/批次', 'UMI count', 'U2OS 24h；四剂量；DMSO', '额外表达监督/结构泛化', '两个RData已下载并实际检查', '逐批导出/公共化合物结构与注释核对'],
        ['Tahoe', 'Tahoe公开Hugging Face发布', 'https://huggingface.co/datasets/tahoebio/Tahoe-100M', 'data/raw/Tahoe_100M', '单细胞表达分片', 'gene token及count；首CLS非基因', '各sample处理；元数据覆盖不等于本地实际覆盖', '先作格式验证，后选择性表达预训练', '1/3388分片28225行＋3元数据；非全库', '补gene token映射并选所需条件分片']
    ]
    # The DRUG-seq URL is taken from the project's actual download provenance below.
    manifest = json.loads((ROOT / 'data/download_manifest.json').read_text())
    drugseq = [r for r in manifest if r['dataset'] == 'Novartis_DRUGseq']
    if drugseq:
        rows[11][2] = '|'.join(r['url'] for r in drugseq)
    pd.DataFrame(rows, columns=columns).to_csv(OUT / 'dataset_registry.csv', index=False)
    claims = [
        ['联合预测响应和功能', 'HMS有同项目条件对应，矩阵已提取', '未训练；性能收益未证明', 'H1独立任务vs共享任务，未见药物', '两主要终点的效应及区间支持相应主张'],
        ['实测早期表达增加后期预测信息', '主队列3h/24h与72h观测', '有观测但信息增益未评价', 'H2给药前vs增加实测3h', '额外测量场景分开报告并有泛化增益'],
        ['预测表达改善功能预测', '可训练中间响应目标', '未实施交叉拟合比较', 'H3串联vs匹配容量直接模型', '控制数据/容量后收益稳定；解释为表示收益'],
        ['共享/专属表示有生物意义', '已有MIX-seq关联分解先例', '当前无本项目模型/机制证据', 'H4统计/神经分解、探针和外部验证', '可重复关联可主张；因果机制另需干预'],
        ['状态加delta构成后续状态', 'L3可支持同平台参考路线', 'CCLE logTPM与L5 z相加缺乏量纲依据', 'H6同平台状态定义与时间预测', '一致尺度和独立时间验证'],
        ['跨细胞泛化', '五核心细胞；扩展数据更多身份待核', '五细胞数据准备；未训练', 'LOCO及严格全局去重/预训练边界', '逐细胞结果支持且说明范围'],
        ['跨平台泛化', '多平台实际数据已下载', '外部矩阵处理和评分未完成', 'H5/H6来源留出和外部表达', '独立实验不重用原曲线且可复现'],
        ['序贯策略/世界模型', '单药时间点和既有相关文献', '目前没有本项目序贯模型与轨迹证据', '后续顺序/洗脱/单药对照及多步验证', '对应轨迹和独立实验后升级主张']
    ]
    pd.DataFrame(claims, columns=['proposal_claim', 'current_evidence', 'actual_status', 'required_comparison', 'claim_upgrade_condition']).to_csv(OUT / 'claims_evidence.csv', index=False)
    issues = [
        ['P0', '正文与NC具体版本', '本地找到初步方案和图，已异步询问正文路径/NC指代', '没有已提供的正式正文或审稿意见', '提供后逐项映射；当前继续初步方案与暂定官方规范设计', '论文结论与合规判断'],
        ['P0', '主队列表达生物QC', 'ID/数值/库板已核；81单实例且对应相关指标缺失', 'L3处理/vehicle分布和可靠性未逐条判定', '形成逐签名QC、对照异常及敏感性队列', '入组和训练发布'],
        ['P1', '跨库板差异', '239重复条件；最大差集中CGP60474低剂量；原孔重算可复现且三批均存在，已排除本项目聚合/比例计算造成', '实际物质/暴露/协议或原注释等原因尚未区分', '保留库板并追原实例对照与化合物批次；核查两源protocol', '聚合和稳健性'],
        ['P1', 'MCF10A背景', '正确ModelID；24Q4无profile/表达；四独立RNA已得', '24Q4上游未覆盖原因官方未说明；培养条件和4基因对应待核', '查稳定ID/别名和原培养协议；独立背景路线与六细胞L3方案比较', '第六细胞入组'],
        ['P1', 'WH-4-025结构', 'PhaseII实际SMILES/InChIKey为restricted', '公开同一物质结构待找到与核验', '追原论文/公共化学登记/作者注释；不凭近名替代', '第107药结构输入'],
        ['P1', 'LINCS剂量实例引用', '扩展审计115签名有一个引用未解析；无已证实剂量冲突', '元数据覆盖/版本或引用问题待定位', '追原inst表/另一发布/实例ID；单独保留已解析与未解析', '扩展剂量对应'],
        ['P1', '扩展身份', '名称对应48模型，47有当前背景；SK-LU-1只有WES profile', '稳定身份及替代RNA来源待补', '核Cellosaurus/供应商/STR与SK-LU-1独立RNA覆盖', '扩展可用规模'],
        ['P1', '药物身份冲突', '5个同BRD不同结构；afatinib同名不同ID/立体声明', '盐/立体/批次/登记原因待核', '保存原结构，逐一查源登记与母体；记录合并或分开证据', '跨库合并与全局划分'],
        ['P1', 'GDSC归一化', '两版原孔和assay已核；3个GDSC1板缺NC0/1', '替代control及正式剂量目标未导出', '核官方归一化和原control字段，复现小批再扩大', '功能长表'],
        ['P1', '外部数据重用', 'MIXseq sens来自GDSC/PRISM', '各原curve与训练重叠尚未逐一排除', '建立原实验来源键，移除重叠或单列重用验证', '独立功能证据'],
        ['P2', 'DRUGseq结构注释', '两RData计数和批次已检查', '公共化合物结构/注释尚未完整接入', '查作者处理注释，逐批导出并核ID', '额外表达训练'],
        ['P2', 'Tahoe格式与范围', '1/3388分片＋3元数据；CLS已识别', 'gene token表和实际条件覆盖待核', '下载映射并核token，再按需要选分片', '预训练与外部表达规模'],
        ['P1', '实测landmark与delta清单', '两期实测978相同；delta表交集785，原因层级已定位为不同集合', '不能将同978行解释成同一基因集合', '所有输出按gene ID与pr_is_lm对齐，注册集合指纹', '避免静默列错配']
    ]
    pd.DataFrame(issues, columns=['priority', 'issue', 'verified_evidence', 'remaining_question', 'next_action', 'affected_result']).to_csv(OUT / 'issue_register.csv', index=False)
    summary = {'original_manifest_files': len(inventory), 'all_files_exist_and_match_recorded_size': True,
               'original_raw_bytes': int(inventory.actual_bytes.sum()), 'original_raw_decimal_GB': float(inventory.actual_bytes.sum()/1e9),
               'raw_directories': int(inventory.dataset.nunique()), 'explanatory_dataset_items': len(rows),
               'additional_MCF10A_RNAseq_files': 4, 'inventory_check': 'Filesystem stat compared with prior download records; no new full-raw checksum scan.',
               'model_fits_this_design_turn': 0}
    (OUT / 'delivery_snapshot.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
