# -*- coding: utf-8 -*-
import gzip, csv
from collections import Counter

samp = Counter(); ct = Counter(); n = 0
with gzip.open(r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\senescence_project\results\fullrun\x23_scores_GSE269748.csv.gz", "rt") as f:
    r = csv.DictReader(f)
    for row in r:
        n += 1
        samp[row["sample"]] += 1
        ct[row["cell_type"]] += 1
print("total cells", n)
print("samples", len(samp))
print(sorted(k + "=" + str(v) for k, v in samp.items()))
print("cell types", len(ct))
print(sorted(k + "=" + str(v) for k, v in ct.items()))
