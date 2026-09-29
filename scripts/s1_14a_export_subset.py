# s1_14a_export_subset.py — 把 7d/6mo 子集从 npz 转为 R 可读格式
# 输入: results/planA_step35_fm_prep/GSE269748/{_subset_7d6mo.npz, genes.txt, _target_cells.csv}
# 输出: results/s1_269748_7d6mo/{counts.mtx, genes.txt, cells.csv}
#   用 MatrixMarket 文本格式 (R 用 Matrix::readMM 直读); 该子集非零仅数千万级, 可接受
import os, time
import numpy as np, pandas as pd
from scipy.sparse import load_npz
from scipy.io import mmwrite

SEN  = r"D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
INPD = os.path.join(SEN, "results/planA_step35_fm_prep/GSE269748")
OUTD = os.path.join(SEN, "results/s1_269748_7d6mo"); os.makedirs(OUTD, exist_ok=True)

t0 = time.time()
npz = os.path.join(INPD, "_subset_7d6mo.npz")
if not os.path.exists(npz):
    raise SystemExit(f"[s1_14a] 缺 {npz} —— 抽取尚未完成")

M = load_npz(npz).tocsc()
print(f"[s1_14a] 子集 {M.shape} 非零={M.nnz:,} ({time.time()-t0:.0f}s)", flush=True)

genes = [l.strip() for l in open(os.path.join(INPD, "genes.txt"))]
tgt   = pd.read_csv(os.path.join(INPD, "_target_cells.csv"))
assert M.shape[0] == len(genes), f"基因数不符 {M.shape[0]} vs {len(genes)}"
assert M.shape[1] == len(tgt),   f"细胞数不符 {M.shape[1]} vs {len(tgt)}"

# 基因名做 safe 处理 (mmwrite 不支持奇怪字符; 这里只确保非空且唯一)
assert len(set(genes)) == len(genes), "基因名有重复"
tgt.to_csv(os.path.join(OUTD, "cells.csv"), index=False)
with open(os.path.join(OUTD, "genes.txt"), "w") as f:
    f.write("\n".join(genes) + "\n")

# 稀疏 mtx (coordinate, 1-based)
mmwrite(os.path.join(OUTD, "counts.mtx"), M, field="integer")
print(f"[s1_14a] DONE → {OUTD}/counts.mtx  ({time.time()-t0:.0f}s)", flush=True)

# 自检
print(f"  细胞: {M.shape[1]} | 基因: {M.shape[0]}", flush=True)
print(f"  逐时点 x 型:\n{tgt.groupby(['time_point','cell_type']).size().to_string()}", flush=True)
