# -*- coding: utf-8 -*-
"""s1_20_L4_translational.py — L4 转化层: 慢性 p16 维持程序的候选干预
三步:
  ① 构建"慢性维持签名": SCENIC persistent regulon 靶基因 + p16+ 表型升高基因
  ② 药物-靶点交叉: 已知 senomorphic/senolytic 药物靶点 vs 签名基因
  ③ 靶点表达验证: 候选药靶在 p16+ 小胶质中是否表达
产物: results/s1_L4_translational/{candidates.csv, LOG.txt}
"""
import os, json
import numpy as np, pandas as pd

SEN  = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\senescence_project"
OUTD = os.path.join(SEN, "results/s1_L4_translational"); os.makedirs(OUTD, exist_ok=True)
L = []
def log(s=""):
    print(s, flush=True); L.append(s)

# ==================== ① 慢性维持签名 ====================
log("=== ① 构建慢性维持签名 ===")

# 来源 A: SCENIC persistent regulons
scenic_f = os.path.join(SEN, "results/microglia_mechanism/scenic/phaseC_10_classification.csv")
scenic = pd.read_csv(scenic_f)
persistent_tfs = scenic[scenic["class"] == "persistent"]["regulon"].tolist()
log(f"Persistent TFs: {persistent_tfs}")

# 来源 B: p16+ 表型升高基因 (from s1_17b2)
pheno_genes_elevated = {
    "Cst7": 23.2,   # +pp in p16+ vs p16-
    "Cd9": 14.0,
    "Itgax": 9.3,
    "Apoe": -3.3,   # not elevated, but key DAM marker
    "Trem2": 5.8,
}
log(f"p16+ 升高基因: {pheno_genes_elevated}")

# 来源 C: 维持信号配体 (from L1 NicheNet)
maintenance_ligands = ["Mif", "Adam17", "Tnf", "Osm", "Il1b"]
log(f"维持信号配体: {maintenance_ligands}")

# 综合签名 = persistent TFs + 升高基因 + 维持配体
signature_genes = list(set(
    persistent_tfs + list(pheno_genes_elevated.keys()) + maintenance_ligands
))
log(f"慢性维持签名: {len(signature_genes)} 基因 → {sorted(signature_genes)}")

# ==================== ② 药物-靶点交叉 ====================
log("\n" + "=" * 76)
log("=== ② 已知 senomorphic/senolytic 药物靶点 vs 签名 ===")

# 已知干预靶点 → 药物映射 (文献来源, 逐条标注)
DRUG_TARGETS = {
    # Senolytics (PI 裁定不进主线, 仅列出供参考)
    "Navitoclax (ABT-263)": {"targets": ["BCL2", "BCL2L1", "BCL2L11"], "class": "senolytic",
                              "caveat": "PI 裁定不进主线"},
    "Dasatinib+Quercetin": {"targets": ["STAT3", "SRC", "EGFR"], "class": "senolytic",
                             "caveat": "PI 裁定不进主线"},
    # Senomorphics (我们的主方向)
    "Ruxolitinib": {"targets": ["JAK1", "JAK2"], "class": "senomorphic",
                     "rationale": "JAK-STAT 通路 = Stat1/Stat2 persistent regulon 的上游"},
    "Tofacitinib": {"targets": ["JAK1", "JAK3"], "class": "senomorphic",
                     "rationale": "同上, JAK 抑制"},
    "Metformin": {"targets": ["AMPK", "mTOR"], "class": "senomorphic",
                   "rationale": "AMPK 激活 → 抑制 SASP 转录"},
    "Rapamycin": {"targets": ["mTOR"], "class": "senomorphic",
                   "rationale": "mTOR 抑制 → 减少 SASP 翻译"},
    # 代谢检查点 (Nature 2026 SLC25A1 轴)
    "SLC25A1 抑制": {"targets": ["SLC25A1"], "class": "senomorphic",
                      "rationale": "Nat 2026 (PMID 42527602): 乙酰-CoA 供给 = SASP 必要条件"},
    # 抗 CD47/SIRPα (免疫清除)
    "Anti-CD47": {"targets": ["CD47"], "class": "immune clearance",
                   "rationale": "促进衰老细胞被吞噬"},
    # BCL-2 家族选择性
    "A-1155463": {"targets": ["BCL2L1"], "class": "senomorphic (BH3 mimetic)",
                   "rationale": "BCL-XL 选择性抑制 → senescent 细胞应激"},
    # MCL-1 (急性程序靶, 非慢性)
    "S63845": {"targets": ["MCL1"], "class": "MCL-1 inhibitor",
                "rationale": "⚠️ 急性 p21+ 小胶质 Mcl1-high → 清除可能有害 (Valdivieso 2026)"},
}

rows = []
for drug, info in DRUG_TARGETS.items():
    targets = info["targets"]
    # 检查靶点是否在签名中
    in_sig = [t for t in targets if t.upper() in [g.upper() for g in signature_genes]]
    # 检查靶点是否在人类同源基因中 (大小写转换)
    rows.append({
        "drug": drug,
        "targets": ", ".join(targets),
        "class": info["class"],
        "rationale": info.get("rationale", info.get("caveat", "")),
        "targets_in_signature": ", ".join(in_sig) if in_sig else "无直接靶点",
        "n_sig_hits": len(in_sig),
    })
drug_df = pd.DataFrame(rows)
log(drug_df[["drug", "class", "targets_in_signature", "n_sig_hits"]].to_string(index=False))

# ==================== ③ 靶点表达验证 (p16+ 小胶质) ====================
log("\n" + "=" * 76)
log("=== ③ 靶点表达验证 (p16+ 小胶质, 来自 s1_17a2 面板) ===")
# 读入 6mo 面板
panel_f = os.path.join(SEN, "results/s1_phenotype/mg_panel_6mo_v2.csv.gz")
if os.path.exists(panel_f):
    panel = pd.read_csv(panel_f)
    p16pos = panel["Cdkn2a"] > 0
    p16neg = panel["Cdkn2a"] == 0
    log(f"p16+ n={p16pos.sum()} | p16− n={p16neg.sum()}")

    # 检查关键靶点基因在面板中的表达
    check_genes = ["Stat1", "Irf7", "Cst7", "Itgax", "Cd9", "Trem2", "Apoe",
                   "Mcl1", "Bcl2", "Bax", "Mif", "Lyz2"]
    # JAK 通路基因 (不在面板中, 用替代读出)
    log(f"\n{'基因':>10s} {'p16+ 检出%':>10s} {'p16− 检出%':>10s} {'方向':>6s}")
    for g in check_genes:
        if g not in panel.columns: continue
        det_pos = 100 * (panel.loc[p16pos, g] > 0).mean()
        det_neg = 100 * (panel.loc[p16neg, g] > 0).mean()
        direction = "↑" if det_pos > det_neg + 2 else ("↓" if det_pos < det_neg - 2 else "=")
        log(f"{g:>10s} {det_pos:10.1f} {det_neg:10.1f} {direction:>6s}")

    # JAK-STAT 通路读出 (Stat1 = JAK 下游)
    stat1_pos = 100 * (panel.loc[p16pos, "Stat1"] > 0).mean()
    stat1_neg = 100 * (panel.loc[p16neg, "Stat1"] > 0).mean()
    log(f"\nJAK-STAT 通路 (Stat1): p16+ {stat1_pos:.1f}% vs p16− {stat1_neg:.1f}%")
    if stat1_pos > stat1_neg + 5:
        log("  → JAK 抑制剂 (ruxolitinib/tofacitinib) 靶向 p16+ 细胞的合理性: 有(Stat1 阳性)")
    else:
        log("  → JAK 抑制剂: Stat1 在 p16+ 无富集, 需要用 persistent regulon 活性而非检出率判")
else:
    log("⚠️ 面板文件不存在, 跳过表达验证")

# ==================== 输出 ====================
drug_df.to_csv(os.path.join(OUTD, "candidates.csv"), index=False)
log(f"\n[L4] DONE → {OUTD}")
with open(os.path.join(OUTD, "LOG.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
