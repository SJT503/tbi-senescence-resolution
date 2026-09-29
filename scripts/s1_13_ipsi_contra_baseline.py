# s1_13_ipsi_contra_baseline.py — 同侧主导的分母校正 (排除"基线低→升幅大"假象)
# 问题: s1_12 报 Oligo 同侧升幅 +615% 排名第1, 但 Oligo 的 naive 基线仅 1.31%(全场最低之一)
#       → 升幅大可能纯粹是分母小, 不是"对侧特别干净"
# 三种校正口径:
#   ① 原始升幅 (ipsi-contra)/contra            [s1_12 已报, 未校正]
#   ② 基线归一化: (ipsi-naive)/naive 与 (contra-naive)/naive, 报比值
#   ③ 绝对差 ipsi-contra (不看分母, 看绝对值)
#   ④ 层级校正: 把该型的 ipsi/contra 与全场中位比值比较 (排除"该型整体偏高/偏低")
# 数据: processed_data/s1_24h_pilot/_covary_export.csv.gz
# 用法: python s1_13_ipsi_contra_baseline.py
# 产物: results/s1_crosscelltype_24h/{ipsi_contra_corrected.csv, LOG_BASELINE.txt}
import os
import numpy as np, pandas as pd

SEN  = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\senescence_project"
INP  = os.path.join(SEN, "results/s1_crosscelltype_24h/permouse_alltypes.csv")
OUTD = os.path.join(SEN, "results/s1_crosscelltype_24h")
L = []
def log(s=""):
    print(s, flush=True); L.append(s)

pm = pd.read_csv(INP)
log(f"[s1_13] 逐鼠逐型表 {len(pm)} 行 | 类型 {pm.cell_type.nunique()}")

PAIRS = [("RNASEQ17L","RNASEQ17R"), ("RNASEQ18L","RNASEQ18R"), ("RNASEQ19L","RNASEQ19R")]
NAIVE = ["RNASEQ04","RNASEQ08","RNASEQ10","RNASEQ16","RNASEQ23","RNASEQ24","RNASEQ25"]

def val(ty, sid):
    r = pm[(pm.cell_type == ty) & (pm["sample"] == sid)]
    return r.det.values[0] if len(r) else np.nan

rows = []
for ty in sorted(pm.cell_type.unique()):
    if ty in ("UNMATCHED", "Ambiguous"):
        continue
    pairs = [(val(ty, a), val(ty, b)) for a, b in PAIRS]
    pairs = [(x, y) for x, y in pairs if not (np.isnan(x) or np.isnan(y))]
    nv = pm[(pm.cell_type == ty) & (pm["sample"].isin(NAIVE))].det.values
    if len(pairs) < 2 or len(nv) < 2:
        continue
    ip = np.array([x for x, _ in pairs]); ct = np.array([y for _, y in pairs])
    naive = float(np.mean(nv))
    eps = 1e-9
    # ① 原始升幅
    lift_raw = np.nanmean([100 * (x - y) / y if y > 0 else np.nan for x, y in pairs])
    # ② 基线归一化: 各自相对 naive 的富集倍数, 再取比
    ip_enr = ip / max(naive, eps)          # 同侧富集倍数
    ct_enr = ct / max(naive, eps)          # 对侧富集倍数
    ratio_norm = float(np.mean(ip_enr / np.maximum(ct_enr, eps)))
    # ③ 绝对差
    diff_abs = float(np.mean(ip - ct))
    # ④ 归一化后仍看"对侧是否被压住": contra/naive 应接近 1
    rows.append(dict(cell_type=ty, n_pairs=len(pairs), naive=naive,
                     ipsi=float(np.mean(ip)), contra=float(np.mean(ct)),
                     lift_raw_pct=lift_raw,
                     ipsi_enrich=float(np.mean(ip_enr)),
                     contra_enrich=float(np.mean(ct_enr)),
                     ratio_norm=ratio_norm,
                     diff_abs=diff_abs,
                     contra_over_naive=float(np.mean(ct_enr))))
d = pd.DataFrame(rows)
d["rank_lift_raw"]   = d.lift_raw_pct.rank(ascending=False)
d["rank_ratio_norm"] = d.ratio_norm.rank(ascending=False)
d["rank_diff_abs"]   = d.diff_abs.rank(ascending=False)
d = d.sort_values("ratio_norm", ascending=False).reset_index(drop=True)
d.to_csv(os.path.join(OUTD, "ipsi_contra_corrected.csv"), index=False)

log("\n" + "=" * 100)
log("同侧主导 — 四种口径对比 (按基线归一化比值排序)")
log(f"{'cell_type':>26s} {'n':>2s} {'naive':>7s} {'ipsi':>8s} {'contra':>8s} "
    f"{'原始升幅%':>10s} {'ipsi富集':>9s} {'contra富集':>10s} {'归一化比':>9s} {'绝对差':>8s}")
for _, r in d.iterrows():
    log(f"{r.cell_type:>26s} {int(r.n_pairs):>2d} {r.naive:7.2f} {r.ipsi:8.2f} {r.contra:8.2f} "
        f"{r.lift_raw_pct:10.1f} {r.ipsi_enrich:9.2f} {r.contra_enrich:10.2f} "
        f"{r.ratio_norm:9.2f} {r.diff_abs:8.2f}")

log("\n" + "=" * 100)
log("排名对比 (第1名是谁):")
for col, name in [("lift_raw_pct","① 原始升幅"), ("ratio_norm","② 基线归一化比"), ("diff_abs","③ 绝对差")]:
    top = d.loc[d[col].idxmax()]
    log(f"  {name:16s}: 第1 = {top.cell_type:26s} (值 {top[col]:.2f})")
    ol = d[d.cell_type == "Oligo"]
    if len(ol):
        rk = int(d[col].rank(ascending=False)[d.cell_type == "Oligo"].values[0])
        log(f"  {'':16s}  Oligo = {ol[col].values[0]:.2f} → 排名 {rk}/{len(d)}")

log("\n" + "=" * 100)
log("关键判读 (Oligo 是否真的是'对侧特别干净'):")
ol = d[d.cell_type == "Oligo"]
if len(ol):
    o = ol.iloc[0]
    log(f"  Oligo naive={o.naive:.2f} | 同侧富集={o.ipsi_enrich:.2f}x | 对侧富集={o.contra_enrich:.2f}x")
    log(f"  → 对侧相对 naive 富集 {o.contra_enrich:.2f}x ", )
    log(f"     {'(对侧被压住, 接近1 → 支持局部化)' if o.contra_enrich < 2 else '(对侧也明显升高 → 非局部化)'}")
    med_ctrl = d[d.cell_type != 'Oligo'].contra_enrich.median()
    log(f"  其他型对侧富集中位 = {med_ctrl:.2f}x")
    log(f"  → 归一化比排名: Oligo 第 {int(d.ratio_norm.rank(ascending=False)[d.cell_type=='Oligo'].values[0])}/{len(d)}")

with open(os.path.join(OUTD, "LOG_BASELINE.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
log(f"\n[s1_13] DONE → {OUTD}")
