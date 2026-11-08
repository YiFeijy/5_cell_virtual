#!/usr/bin/env python3
"""Register selected primary works and claim overlap; not a systematic review."""
from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / '项目总设计_20261106/相关工作与创新性复核_20261106'


def main():
    OUT.mkdir(exist_ok=True)
    fields = ['id', 'work', 'year', 'venue_and_status', 'doi_or_identifier', 'primary_url', 'input_and_supervision', 'output_or_question', 'overlap_with_current_design', 'important_difference', 'verification_depth', 'evidence_locator']
    rows = [
        [1, 'Niepel et al.', 2017, 'Nature Communications;期刊论文', '10.1038/s41467-017-01383-w', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC5662764/', '同项目L1000和初始/72h平行计数', '表达与生长/存活的关系和分歧', '当前主数据及科学现象的直接先例', '原研究不是当前跨药物联合网络', '原文摘要/Introduction/Results；本地数据另行核过', 'Abstract; Results: High-dimensional drug response profiling; transcriptional/phenotypic response classes'],
        [2, 'Szalai et al.', 2019, 'Nucleic Acids Research;期刊论文', '10.1093/nar/gkz805', 'https://doi.org/10.1093/nar/gkz805', 'LINCS实测签名和近剂量CTRP等', '由表达预测活性；通用死亡/增殖签名及混杂', '表达到功能与共同响应信号', '功能关系及统计模型不等于当前基线到双目标', '原文摘要/方法检索内容；已有本地作者notebook核查', '作者cell_death_1_preprocess.ipynb; Methods linear ridge regression'],
        [3, 'Dr.VAE', 2019, 'Bioinformatics;期刊论文', '10.1093/bioinformatics/btz158', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC6761940/', '处理前/6h药后表达与CTRP敏感性，允许缺标签', '潜在扰动映射及药敏分类', '响应辅助功能预测和半监督状态建模', '每药单独模型；AAC二分类，不是跨药结构/连续剂量头', '原文Methods/Results', 'Methods generative model; Results 3.2 and 3.4; each of 26 drugs'],
        [4, 'MIX-seq', 2020, 'Nature Communications;期刊论文', '10.1038/s41467-020-17440-w', 'https://doi.org/10.1038/s41467-020-17440-w', 'scRNA药物响应与GDSC/PRISM AUC等', '活性相关/无关表达；早期表达预测药敏', '分解故事与早期表达价值的直接先例', '统计斜率/截距；曲线敏感性，不是独立同剂量72h计数', '原文Fig.2/3和方法检索内容；本地README/代码记录', 'Fig.2 viability-related/-independent; Fig.3 prediction; Drug sensitivity data'],
        [5, 'MultiDCP', 2022, 'PLOS Computational Biology;期刊论文', '10.1371/journal.pcbi.1010367', 'https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010367', '基线表达/结构/剂量；异构表达及拟合曲线标签', '扰动表达与剂量相关活性联合预测', '最直接的双任务先例', 'PharmacoDB曲线估计；当前HMS是平行条件实际计数', '原文Architecture及功能数据方法', 'Architecture of multi-task dose-dependent chemical phenomics model; Cell viability prediction dataset'],
        [6, 'scTherapy', 2024, 'Nature Communications;期刊论文', '10.1038/s41467-024-52980-5', 'https://doi.org/10.1038/s41467-024-52980-5', 'LINCS/药物特征和PharmacoDB曲线插值', '功能预测与患者药物优先排序', '表达功能对接/目标浓度估计', '不是同信息条件下联合生成表达与计数', '原文方法检索片段及既有项目核查', 'LINCS 2020 concentrations matched with PharmacoDB; LightGBM outcome interpolated cell viability'],
        [7, 'SequenTx', 2026, 'Nature Machine Intelligence;期刊论文', '10.1038/s42256-026-01192-1', 'https://doi.org/10.1038/s42256-026-01192-1', 'LINCS表达、CTRP功能及序贯实验', '状态/活性模型与强化学习序贯策略', '长期转移/功能/序贯设想的直接先例', '当前仅单药数据准备，没有相应序贯结果', '原文摘要/数据代码说明；作者预处理及已有本地代码', 'Code availability; DELTA-TJ/SequenTx/reproducity/CellViabilityModel/01_data_preprocess.py'],
        [8, 'Domain Separation Networks', 2016, 'NeurIPS;会议论文', 'NeurIPS2016:45fbc6d3e05ebd93369ce542e8f2322d', 'https://proceedings.neurips.cc/paper_files/paper/2016/hash/45fbc6d3e05ebd93369ce542e8f2322d-Abstract.html', '源域/目标域图像', '共享和各域私有表示及重建', '共享/专属结构的一般先例', '图像域适配，非表达/活性任务', '原始会议摘要及原论文索引', 'Abstract; private/shared subspaces and reconstruction'],
        [9, 'CODE-AE', 2022, 'Nature Machine Intelligence;期刊论文', '10.1038/s42256-022-00541-0', 'https://doi.org/10.1038/s42256-022-00541-0', '细胞系/患者表达，细胞系药敏标签', '共享/私有编码、对齐和临床药敏迁移', '生物药敏中的共享/私有结构', '分离数据域；本方案拟分任务', '原文Overview/Fig.1/Methods检索正文', 'Fig.1 shared/private encoders; soft orthogonality; fine-tuning binary AUC labels'],
        [10, 'CPA', 2023, 'Molecular Systems Biology;期刊论文', '10.15252/msb.202211517', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC10258562/', '单细胞表达和扰动/剂量/协变量', '未见条件的转录响应生成', '背景/干预分离、剂量时间条件化', '主要目标表达，非独立计数双任务', '原文摘要/Introduction及方法概述', 'Abstract; modular chemical representation; unseen dosage/time/cell types'],
        [11, 'biolord', 2024, 'Nature Biotechnology;期刊论文', '10.1038/s41587-023-02079-x', 'https://doi.org/10.1038/s41587-023-02079-x', '单细胞多组学及已知属性', '已知/未知属性分解与状态生成', '分解与跨状态预测先例', '属性分解非功能相关/无关两任务', '原始摘要及简述；未完整复现方法', 'Abstract: known and unknown attributes; unseen drugs/genetic perturbations'],
        [12, 'scDisInFact', 2024, 'Nature Communications;期刊论文', '10.1038/s41467-024-45227-w', 'https://doi.org/10.1038/s41467-024-45227-w', '多批次多条件scRNA表达', '共享/非共享表示、整合与条件预测', '单细胞共享/特异结构先例', '不等于表达与剂量功能联合预测', '原文摘要/方法索引片段；未逐项读实现', 'Shared/unshared encoders and condition-associated genes'],
        [13, 'TranSiGen', 2024, 'Nature Communications;期刊论文', '10.1038/s41467-024-49620-3', 'https://doi.org/10.1038/s41467-024-49620-3', '同板DMSO/处理表达与化合物表示', '基线/处理表达分布和映射', '药物条件表达预训练与状态映射', '同平台参考与跨研究CCLE背景不同', '原文Overview和paired preprocessing方法检索内容', 'Overview of TranSiGen; same-plate DMSO paired X1-X2'],
        [14, 'XPert', 2026, 'Nature Machine Intelligence;期刊论文', '10.1038/s42256-025-01165-w', 'https://doi.org/10.1038/s42256-025-01165-w', '基线表达、化学/知识表示、剂量时间', '药后表达和变化量', '条件交互、双分支和剂量时间建模', '两个输出均是转录目标，不是独立功能头', '原文Overview/Fig.1/Methods检索正文', 'Base self-attention vs perturb cross-attention; mdmt; x_pert and x_deg'],
        [15, 'Cell-JEPA', 2026, 'arXiv v1预印本；本轮未确认其他发表版本', 'arXiv:2602.02093v1', 'https://arxiv.org/html/2602.02093v1', '掩码/完整scRNA；遗传扰动下游数据', '细胞表示和药后绝对状态/变化量评价', '细胞JEPA已有先例', '已读扰动实验是K562 CRISPR；非药物与独立功能联合；delta未获益', '原文HTML方法、结果及表3/4', 'Sections 2.2/2.3/3.4; Table4 delta metrics'],
        [16, 'MultiVCDiff', 2026, 'Cell Reports Methods;期刊论文', '10.1016/j.crmeth.2026.101459', 'https://pubmed.ncbi.nlm.nih.gov/42167225/', '化学或遗传扰动；表达/形态训练语料', '表达与形态联合生成', '多模态联合输出已有方法', '形态不是独立剂量功能测量', '原始摘要及Fig.1方法图注', 'PubMed Abstract and Figure1 paired latent encoders/conditional denoising'],
        [17, 'BRiCE', 2026, 'Scientific Reports;期刊论文', '10.1038/s41598-026-54526-9', 'https://doi.org/10.1038/s41598-026-54526-9', 'T/B细胞表达与克隆扩增', '响应细胞识别和强度基准', '表达与丰度区别及验证问题', '免疫克隆扩增，不是药物计数预测主任务', '原文摘要、Introduction及验证设计', 'RC and Res; clonotype expansion as validation'],
        [18, 'Counting cells can accurately predict small-molecule bioactivity benchmarks', 2026, 'Nature Communications;期刊论文', '10.1038/s41467-026-68725-5', 'https://doi.org/10.1038/s41467-026-68725-5', '细胞数与小分子生物活性基准', '简单计数基线与表型基准混杂', '必须区分一般数量/毒性与复杂表示收益', '未据未读全文扩写其数据与实现', '仅原期刊检索摘要片段；全文跳转未成功', 'Original journal abstract indexed excerpt']
    ]
    assert len(rows) == 18
    pd.DataFrame(rows, columns=fields).to_csv(OUT / 'paper_comparison.csv', index=False)
    claims = [
        ['表达与剂量功能联合预测', '已有直接先例', 'MultiDCP;Dr.VAE', '统一数据预算/标签与最近模型直接比较，提出额外可验证问题'],
        ['共享与专属网络', '结构已有；当前功能分工尚未识别', 'DSN;CODE-AE;scDisInFact', '消融/探针/稳定性/外部验证；排查共享支路被绕过'],
        ['表达变化与生长功能存在分歧', '当前主数据原论文已有', 'Niepel2017;Szalai2019;MIXseq2020', '超出重复现象的新规律、预测或解释证据'],
        ['早期表达预测后期功能', '已有直接先例', 'MIXseq2020;Dr.VAE', '严格信息预算、时间成本及未见药物/来源的增量价值'],
        ['剂量时间和潜在转移', '已有直接相关工作', 'CPA;TranSiGen;XPert;SequenTx', '一致尺度与实际时点/轨迹；新映射应有特定检验'],
        ['JEPA用于细胞', '已有预印本先例；不能泛称首次', 'Cell-JEPA', '检查干预条件化是否提升变化量与功能，区别背景保持'],
        ['细胞世界模型/序贯策略', '已有相关工作；当前项目未验证', 'SequenTx', '相应序贯轨迹与独立验证后再升级主张'],
        ['同项目真实计数及可追溯对接', '有用验证条件，非自动新方法', 'Niepel2017;本项目数据审计', '证明对可复现结论、可靠性与模型比较的实际贡献']
    ]
    pd.DataFrame(claims, columns=['proposed_claim', 'assessment', 'closest_works', 'required_additional_evidence']).to_csv(OUT / 'novelty_claim_audit.csv', index=False)
    log = {
        'date': '2026-11-06', 'scope': 'Targeted related-work and novelty audit of the current proposed architecture; selected primary papers, not an exhaustive systematic review.',
        'selected_works': len(rows), 'queries_and_axes': [
            'joint prediction drug induced gene expression and dose dependent cell viability multi task learning',
            'shared private disentangled representation gene expression drug response CODE-AE',
            'joint embedding predictive architecture JEPA cellular perturbation gene expression drug',
            'MIX-seq viability dependent independent transcriptional responses early predict sensitivity',
            'MultiDCP joint 2025 / 2026', 'Exact titles and DOIs; original paper methods/figures and author code'],
        'sources': ['Original journal pages/indexed primary text', 'PMC/PubMed', 'arXiv original manuscript', 'NeurIPS proceedings', 'Author GitHub'],
        'connector_failure': {'service': 'Elicit search_papers', 'initial_semantic_requests': 3, 'keyword_retry_requests': 1, 'returned_error': 'INVALID_ARGUMENT', 'usable_results': 0, 'fallback': 'Web search and original primary sources'},
        'access_limits': ['Some Nature/PMC direct opens redirected to login/captcha; indexed original text and existing local evidence used with verification depth recorded.', 'For selected newer boundary papers only abstract/figure or indexed excerpt checked; no full-method claims.', 'Cell-JEPA record checked as arXiv v1; no claim that it has no other publication.'],
        'no_new_model_training': True, 'no_author_benchmark_reproduction': True,
        'local_prior_evidence': ['数据可行性调研_20261105/literature_review.md', '数据可行性调研_20261105/evidence/cell_death_1_preprocess.ipynb', '数据可行性调研_20261105/evidence/SequenTx_viability_preprocess.py', 'data/raw/MIXseq/README.txt'],
        'files': ['文献梳理与方案判断.md', 'paper_comparison.csv', 'novelty_claim_audit.csv']}
    (OUT / 'search_and_evidence_log.json').write_text(json.dumps(log, ensure_ascii=False, indent=2) + '\n')
    print(f'Registered {len(rows)} selected works and {len(claims)} claim assessments.')


if __name__ == '__main__':
    main()
