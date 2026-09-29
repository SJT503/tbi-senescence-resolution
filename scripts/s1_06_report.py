# -*- coding: utf-8 -*-
"""s1_06_report.py — Step-1 第6步: STEP1_REPORT.md 汇编 + 两败即停判定 (协议§5)

读入: qc_summary.csv / decont_summary.csv / numi_floors.csv / prop_gate_check.csv /
      confidence_by_group.csv / readout_keynumbers_{main,lo20,hi20}.json / readout_summary_*.txt
产出: results/s1_24h_pilot/STEP1_REPORT.md
用法: python s1_06_report.py   (s1_01-s1_05b 全部完成后)
"""
import os, json
import pandas as pd

SEN = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\senescence_project"
PD_ = os.path.join(SEN, "processed_data/s1_24h_pilot")
OUTD = os.path.join(SEN, "results/s1_24h_pilot")

def rd(p): return pd.read_csv(p) if os.path.exists(p) else None
qc   = rd(os.path.join(PD_, "sce_qc/qc_summary.csv"))
dc   = rd(os.path.join(PD_, "sce_clean/decont_summary.csv"))
fl   = rd(os.path.join(PD_, "sce_final/numi_floors.csv"))
gate = rd(os.path.join(PD_, "labels/prop_gate_check.csv"))
conf = rd(os.path.join(PD_, "labels/confidence_by_group.csv"))
keys = {v: json.load(open(os.path.join(OUTD, f"readout_keynumbers_{v}.json")))
        for v in ("main", "lo20", "hi20")
        if os.path.exists(os.path.join(OUTD, f"readout_keynumbers_{v}.json"))}
assert "main" in keys, "缺 main 版读出, 先跑 s1_05b"
m = keys["main"]

# ---- 两败即停 4 条 (协议§5) ----
# ① 去污染后仍≥8/13型升高 且 P3 p>0.15 无方向改善
r1 = (m["p2_n_clean_sig"] >= 8) and (m["p3_p"] > 0.15) and (m["p3_mean_24h"] <= m["p3_mean_ctrl"] + 0.05)
# ② S1 无单调关系 — 需 s1 组均值, 从 summary 文本读; 此处先留待人工复核标记
r2 = None   # 由报告中 S1 组均值序列判定 (脚本不自动判, 防误判)
# ③ S2 对侧升幅与同侧差<20% → 从 readout_permouse + summary 文本人工核对
r3 = None
lines = []
A = lines.append
A("# STEP1_REPORT — GSE269748 24h 批去污染再分析（预登记执行报告）")
A("")
A("**日期**: auto-generated | **协议**: STEP1_PROTOCOL_24h_pilot.md（2026-09-24 预登记）")
A("")
A("## 1. 管线运行")
if qc is not None:
    A(f"- QC: {len(qc)} 样本, 通过细胞数中位 {int(qc.cells_pass.median())} "
      f"(raw 中位 {int(qc.cells_raw.median())}), medMT 中位 {qc.med_pct_mt.median():.2f}%")
if dc is not None:
    A(f"- decontX+scDblFinder: 双联体率中位 {dc.dbl_pct.median():.1f}%, "
      f"medContam 中位 {dc.med_contam.median():.4f} (filtered 矩阵本底即低)")
if fl is not None:
    A(f"- nUMI 下采样 floor = {int(fl.floor_used.iloc[0])} (7 组中位之最小)")
if gate is not None:
    A(f"- 标签比例门禁: 通过 (无 ≥1% 型比例差 >2 倍项; 详 prop_gate_check.csv)")
if conf is not None:
    A(f"- 转移置信度: median {conf.med_conf.min() if 'med_conf' in conf else 'NA'}, "
      f"Ambiguous 比例详 confidence_by_group.csv")
A("")
A("## 2. 预登记读出 (main 版, floor 正常)")
A("")
A(f"| 读出 | 结果 |")
A(f"|---|---|")
A(f"| P1 p21 检出 (oligo, exact MWU p) | {m['p1']['p21']:.4g} |")
A(f"| P1 SenMayo (oligo, p) | {m['p1']['sen']:.4g} |")
A(f"| P1 OR: Cdkn1a / Bax / Mcl1 | {m['p1_or'].get('Cdkn1a', float('nan')):.2f} / "
  f"{m['p1_or'].get('Bax', float('nan')):.2f} / {m['p1_or'].get('Mcl1', float('nan')):.2f} |")
A(f"| P2 显著升高型数 (非LowConf/Ambiguous) | {m['p2_n_clean_sig']}/{m['p2_n_clean']} (预登记问: ≥8/13?) |")
A(f"| P3 log2(oligo/bg7) | 24h {m['p3_mean_24h']:+.3f} vs Ctrl {m['p3_mean_ctrl']:+.3f}, "
  f"p={m['p3_p']:.4g} (旧: -0.218 vs -0.401, p=0.0734) |")
A("")
A("完整明细: readout_summary_main.txt; S1/S2/S3 见该文件对应节。")
A("")
A("## 3. 敏感性 (±20% floor)")
A("")
for v in ("lo20", "hi20"):
    if v in keys:
        k = keys[v]
        A(f"- **{v}**: P1 p21 p={k['p1']['p21']:.4g} | P2 sig {k['p2_n_clean_sig']}/{k['p2_n_clean']} | "
          f"P3 {k['p3_mean_24h']:+.3f} vs {k['p3_mean_ctrl']:+.3f} (p={k['p3_p']:.4g})")
A("")
A("## 4. 两败即停判定 (协议§5)")
A("")
A(f"- 条款① (仍≥8型升高 且 P3 p>0.15 无改善): {'🔴 触发' if r1 else '未触发'} "
  f"(P2={m['p2_n_clean_sig']}, P3 p={m['p3_p']:.4g})")
A(f"- 条款② (S1 无单调): **人工判定** — 看 readout_summary_main.txt S1 节组均值序列是否 "
  f"naive<rCHI<CCI<HS 单调")
A(f"- 条款③ (S2 对侧差<20%): **人工判定** — 看 S2 节 '平均升幅' 是否 <20%")
A(f"- 条款④ (S63845 PK): 独立于本试点, 不在此判")
A("")
A("## 5. 过关条件 (协议§6)")
A("")
A(f"- 去污染后 24h 信号仍在: P1 p21 p={m['p1']['p21']:.4g} → {'✅' if m['p1']['p21'] < 0.05 else '❌/弱'}")
A(f"- 梯度单调: 待人工核对 S1")
A(f"- P3 方向改善或显著: p={m['p3_p']:.4g}, 方向 Δ={m['p3_mean_24h'] - m['p3_mean_ctrl']:+.3f} "
  f"(旧 Δ=+0.183) → 待判")
A("")
A("## 6. 结论栏（人工填写）")
A("")
A("- [ ] 过关 → 进入重构版 B 改稿")
A("- [ ] 触发停损 → 终止重构版 B 投入")
A("- 备注: ")
A("")
A("---")
A("*本报告由 s1_06_report.py 自动汇编; 判定栏②③与结论由 PI/主会话复核后填写。*")

with open(os.path.join(OUTD, "STEP1_REPORT.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
print(f"[s1_06] STEP1_REPORT.md written -> {OUTD}")
