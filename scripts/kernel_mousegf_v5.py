#!/usr/bin/env python3
"""Mouse-GF ISP: 6mo MG — download Mouse-Geneformer from HuggingFace + run ISP"""
import numpy as np, pandas as pd, pickle, torch, os, warnings, shutil, time, json, urllib.request, glob as globlib
import re, subprocess
warnings.filterwarnings("ignore")

DEV = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Device: {DEV}")
if DEV == "cuda": print(f"GPU: {torch.cuda.get_device_name(0)}")
t0 = time.time()

# === Step 1: Download Mouse-GF from HuggingFace ===
print("\n=== Downloading Mouse-Geneformer from HuggingFace ===")
os.environ["HF_HOME"] = "/kaggle/working/hf"
try:
    from huggingface_hub import hf_hub_download, list_repo_files

    # Try model repo first
    for repo_id in ["MPRG/Mouse-Geneformer"]:
        try:
            files = list_repo_files(repo_id, repo_type="model")
            print(f"Model repo {repo_id} files: {files}")
            for f in files:
                if f.endswith(('.pkl', '.safetensors', '.bin', '.json')):
                    path = hf_hub_download(repo_id=repo_id, filename=f, repo_type="model")
                    print(f"  Downloaded: {f} -> {path}")
        except Exception as e:
            print(f"  Model repo {repo_id}: {e}")

    # Try dataset repo for dictionaries
    for repo_id in ["MPRG/Mouse-Genecorpus-20M"]:
        try:
            files = list_repo_files(repo_id, repo_type="dataset")
            print(f"\nDataset repo {repo_id} files: {files[:20]}")  # first 20
            for f in files:
                # ONLY download exact dictionary/model files, NEVER .arrow corpus files
                if f in ['MLM-re_token_dictionary_v1.pkl', 'mouse_gene_median_dictionary.pkl',
                          'model.safetensors', 'pytorch_model.bin', 'config.json',
                          'MLM_model.safetensors', 'model.bin']:
                    path = hf_hub_download(repo_id=repo_id, filename=f, repo_type="dataset")
                    print(f"  Downloaded: {f} -> {path}")
        except Exception as e:
            print(f"  Dataset repo {repo_id}: {e}")
except ImportError:
    print("huggingface_hub not available, trying direct URL download...")
    # Fallback: direct URL download
    BASE = "https://huggingface.co/datasets/MPRG/Mouse-Genecorpus-20M/resolve/main/"
    for fname in ["MLM-re_token_dictionary_v1.pkl", "mouse_gene_median_dictionary.pkl"]:
        try:
            urllib.request.urlretrieve(BASE + fname, f"/kaggle/working/{fname}")
            print(f"  Downloaded: {fname}")
        except Exception as e:
            print(f"  Failed {fname}: {e}")

# List what we got
print("\n=== Downloaded files ===")
for p in sorted(globlib.glob("/kaggle/working/hf/**/*.pkl", recursive=True) +
                globlib.glob("/kaggle/working/hf/**/*.safetensors", recursive=True) +
                globlib.glob("/kaggle/working/*.pkl", recursive=True)):
    print(f"  {p} ({os.path.getsize(p)} bytes)")

# === Step 2: Load mouse data from GitHub ===
print("\n=== Loading mouse data ===")
DATA_URL = "https://raw.githubusercontent.com/SJT503/tbi-senescence-resolution/main/mg_6mo_mouse_native.npz"
DATA_PATH = "/kaggle/working/mg_6mo_mouse_native.npz"
if not os.path.exists(DATA_PATH):
    urllib.request.urlretrieve(DATA_URL, DATA_PATH)
D = np.load(DATA_PATH, allow_pickle=True)
X_raw, genes, p16_pos = D["X"], D["genes"], D["p16_pos"]
print(f"Data: {X_raw.shape[0]} cells × {X_raw.shape[1]} genes | p16+: {p16_pos.sum()}")

# === Step 3: Find and load dictionaries ===
print("\n=== Loading Mouse-GF dictionaries ===")
tokdict = None; median = None

# Search for dictionary files
pkl_files = globlib.glob("/kaggle/working/hf/**/*.pkl", recursive=True) + globlib.glob("/kaggle/working/*.pkl")
for p in pkl_files:
    try:
        d = pickle.load(open(p, "rb"))
        if isinstance(d, dict):
            # Check if it's a token dict (keys are gene names/IDs, values are ints)
            # or a median dict (keys are gene names/IDs, values are floats)
            sample_vals = list(d.values())[:5]
            if all(isinstance(v, int) for v in sample_vals):
                tokdict = d; print(f"  Token dict: {p} ({len(d)} entries)")
            elif all(isinstance(v, (int, float)) for v in sample_vals):
                median = d; print(f"  Median dict: {p} ({len(d)} entries)")
    except: pass

if tokdict is None or median is None:
    # Try to construct from what's available
    print("Could not identify dictionaries from downloaded files.")
    print("All downloaded .pkl files:")
    for p in pkl_files:
        print(f"  {p}")
    # List all files in HF cache
    print("\nAll files in HF cache:")
    for p in sorted(globlib.glob("/kaggle/working/hf/**/*", recursive=True)):
        if os.path.isfile(p): print(f"  {p} ({os.path.getsize(p)} bytes)")
    raise ValueError("Could not load Mouse-GF dictionaries")

# === Step 4: Map mouse genes (symbols -> Ensembl if needed) ===
gene_list = [str(g) for g in genes]
sample_keys = list(tokdict.keys())[:5]
print(f"\nToken dict key format: {sample_keys}")
print(f"Median dict key format: {list(median.keys())[:5]}")

# Majority vote on key format (robust to a few odd entries)
ens_mode = sum(1 for k in tokdict if str(k).startswith("ENSMUSG")) > len(tokdict) // 2
sym2ids = {}
if ens_mode:
    print("Token keys are Ensembl IDs — downloading symbol->Ensembl mapping from GitHub")
    MAP_URL = "https://raw.githubusercontent.com/SJT503/tbi-senescence-resolution/main/mouse_symbol_to_ensembl.tsv"
    MAP_PATH = "/kaggle/working/mouse_symbol_to_ensembl.tsv"
    if not os.path.exists(MAP_PATH):
        urllib.request.urlretrieve(MAP_URL, MAP_PATH)
    for line in open(MAP_PATH):
        p = line.rstrip("\r\n").split("\t")
        if len(p) >= 2:
            sym2ids.setdefault(p[0], []).append(p[1])
    print(f"Mapping table: {len(sym2ids)} unique symbols")
    def to_ens(g):
        # symbol -> Ensembl ID; when a symbol maps to several IDs, prefer the one in tokdict
        ids = sym2ids.get(g)
        if not ids:
            return g  # unmapped stays as-is, will simply not match
        for i in ids:
            if i in tokdict:
                return i
        return ids[0]
    before = sum(1 for g in gene_list if g in tokdict)
    gene_list = [to_ens(g) for g in gene_list]
    after = sum(1 for g in gene_list if g in tokdict)
    print(f"Symbol->Ensembl conversion: tokdict matches {before} -> {after}")
else:
    print("Token keys appear to be gene symbols — using directly")

mapped = sum(1 for g in gene_list if g in tokdict)
print(f"Matches: {mapped}/{len(gene_list)}")

if mapped < 100:
    print("Too few matches — dumping diagnostic info")
    print("First 20 gene_list:", gene_list[:20])
    print("First 20 tokdict keys:", list(tokdict.keys())[:20])
    raise ValueError(f"Only {mapped} genes matched to token dictionary")

# Filter to mappable genes
have = np.array([g in tokdict and g in median for g in gene_list])
X = X_raw[:, have]
gene_use = [g for g, h in zip(gene_list, have) if h]
print(f"Mappable: {X.shape[1]} / {X_raw.shape[1]} genes")

# === Step 5: Load model (weights hosted on Google Drive per official repo README) ===
print("\n=== Loading Mouse-GF model ===")
# Official weights (machine-perception-robotics-group/Mouse-Geneformer README):
#   base model  = Drive id 1gM3gcc3DlNGt5bAcqHbeRxtdMktGeDEg  (6L/4H/256, same arch as human GF base)
#   large model = Drive id 1xKMyFA4JJeRigcJPsU2XNyxEW25Q247u  (12L, E20)
WEIGHT_URL = "https://drive.google.com/uc?id=1gM3gcc3DlNGt5bAcqHbeRxtdMktGeDEg"
WEIGHT_PATH = "/kaggle/working/mouse_geneformer_base.pt"
if not os.path.exists(WEIGHT_PATH):
    subprocess.run(["pip", "install", "-q", "gdown"], check=True)
    import gdown
    gdown.download(url=WEIGHT_URL, output=WEIGHT_PATH, quiet=False)
print(f"Weights file: {WEIGHT_PATH} ({os.path.getsize(WEIGHT_PATH)} bytes)")

def load_state_any(path):
    """Load weights from .pt/.bin/.pth/.safetensors, or from an archive zip that bundles them."""
    import zipfile
    if zipfile.is_zipfile(path):
        names = zipfile.ZipFile(path).namelist()
        print(f"Zip detected: {len(names)} entries, first 10: {names[:10]}")
        if "data.pkl" not in names:  # a plain archive, not a torch-serialized file
            exdir = "/kaggle/working/mousegf_extract"
            zipfile.ZipFile(path).extractall(exdir)
            cands = [p for p in globlib.glob(exdir + "/**/*", recursive=True)
                     if os.path.splitext(p)[1] in (".pt", ".bin", ".pth", ".safetensors")
                     and os.path.getsize(p) > 1_000_000]
            print(f"Weight candidates inside archive: {cands}")
            if not cands:
                raise ValueError("Zip archive contained no weight file (.pt/.bin/.pth/.safetensors)")
            path = max(cands, key=os.path.getsize)
    if path.endswith(".safetensors"):
        from safetensors.torch import load_file
        return load_file(path), path
    st = torch.load(path, map_location="cpu", weights_only=False)
    return st, path

state, WEIGHT_PATH = load_state_any(WEIGHT_PATH)
if hasattr(state, "state_dict"):
    state = state.state_dict()
if not isinstance(state, dict):
    raise ValueError(f"Unexpected weight format: {type(state)}")
print(f"Checkpoint keys: {len(state)}, first 5: {list(state.keys())[:5]}")

# Normalize prefix (checkpoint may be saved from BertModel without 'bert.' prefix)
if not any(k.startswith("bert.") for k in state) and not any(k.startswith("cls.") for k in state):
    state = {("bert." + k): v for k, v in state.items()}

# Infer architecture from the checkpoint itself (no guessing)
emb_key = next(k for k in state if k.endswith("word_embeddings.weight"))
vocab_size, hidden_size = state[emb_key].shape[0], state[emb_key].shape[1]
layer_ids = {int(m.group(1)) for k in state if (m := re.match(r"bert\.encoder\.layer\.(\d+)\.", k))}
n_layers = (max(layer_ids) + 1) if layer_ids else 6
inter_key = next((k for k in state if k.endswith("intermediate.dense.weight")), None)
inter_size = state[inter_key].shape[0] if inter_key is not None else hidden_size * 2
print(f"Inferred arch: vocab={vocab_size}, hidden={hidden_size}, layers={n_layers}, intermediate={inter_size}")

from transformers import BertForMaskedLM, BertConfig
config = BertConfig(vocab_size=vocab_size, hidden_size=hidden_size,
                    num_hidden_layers=n_layers, num_attention_heads=4,
                    intermediate_size=inter_size, max_position_embeddings=2048,
                    hidden_act="silu")
model = BertForMaskedLM(config)
try:
    model.load_state_dict(state, strict=True)
    print("State dict loaded (strict)")
except RuntimeError as e:
    model.load_state_dict(state, strict=False)
    print(f"State dict loaded (non-strict): {e}")

# Guard against silently running a randomly-initialized network
ref = model.state_dict()
matched = sum(1 for k in ref if k in state and hasattr(state[k], "shape") and state[k].shape == ref[k].shape)
print(f"Parameter tensors matched: {matched}/{len(ref)}")
if matched < 0.8 * len(ref):
    raise ValueError(f"Only {matched}/{len(ref)} tensors matched — checkpoint/model mismatch, aborting")

model.eval().to(DEV)
n_params = sum(p.numel() for p in model.parameters())
print(f"Model loaded: {n_params/1e6:.1f}M params on {DEV}")

# === Step 6: ISP (same logic as human GF version) ===
TARGETS_MOUSE = ['Cdkn2a', 'Cdkn1a', 'Mif', 'Adam17', 'Stat1', 'Stat2', 'Irf7',
                 'Il1b', 'Osm', 'Cst7', 'Itgax', 'Mcl1', 'Ppia', 'B2m', 'Actb']

def tokenize_batch(X_batch, gene_batch, max_len=2048):
    seqs = []
    for c in range(X_batch.shape[0]):
        vals = [X_batch[c, j] / max(median.get(gene_batch[j], 1), 1e-9)
                for j in range(len(gene_batch))]
        order = np.argsort(vals)[::-1][:max_len]
        seqs.append([tokdict[gene_batch[j]] for j in order if gene_batch[j] in tokdict])
    return seqs

MAX_CELLS = 4000
rng = np.random.RandomState(42)
p16_idx = np.where(p16_pos)[0]
p16_neg_idx = np.where(~p16_pos)[0]
n_neg = min(MAX_CELLS - len(p16_idx), len(p16_neg_idx))
use_idx = np.sort(np.concatenate([p16_idx, rng.choice(p16_neg_idx, n_neg, replace=False)]))
X_use = X[use_idx]; p16_use = p16_pos[use_idx]
print(f"\nUsing {len(use_idx)} cells ({p16_use.sum()} p16+)")

BATCH = 500; all_seqs = []
for b in range(0, len(use_idx), BATCH):
    all_seqs.extend(tokenize_batch(X_use[b:b+BATCH], gene_use))
print(f"Tokenization done ({time.time()-t0:.0f}s)")

def get_emb(model, tokens, dev):
    ids = torch.tensor([tokens], dtype=torch.long).to(dev)
    with torch.no_grad():
        out = model(input_ids=ids, output_hidden_states=True)
    return out.hidden_states[-1][0, 0, :].cpu().numpy()

print("Computing original embeddings...")
orig_embs = []
for i, seq in enumerate(all_seqs):
    orig_embs.append(get_emb(model, seq, DEV))
orig_embs = np.array(orig_embs)
print(f"Embeddings: {orig_embs.shape} ({time.time()-t0:.0f}s)")

p16neg_cent = orig_embs[~p16_use].mean(axis=0)
orig_dist_pos = np.mean([np.linalg.norm(orig_embs[i] - p16neg_cent) for i in range(len(all_seqs)) if p16_use[i]])
print(f"p16+ dist to p16- centroid: {orig_dist_pos:.4f}")

# KO analysis (targets are MGI symbols — convert to Ensembl in ens_mode)
tgt_tokens = {}
for t in TARGETS_MOUSE:
    key = to_ens(t) if ens_mode else t
    tok = tokdict.get(key)
    if tok is not None:
        tgt_tokens[t] = tok
    elif key in gene_use:
        tgt_tokens[t] = tokdict[key]
print(f"\nKO targets: {len(tgt_tokens)}/{len(TARGETS_MOUSE)}")
if len(tgt_tokens) < len(TARGETS_MOUSE):
    missing = [t for t in TARGETS_MOUSE if t not in tgt_tokens]
    print(f"  Missing targets (no token): {missing}")

results = []
for gene, ko_token in tgt_tokens.items():
    print(f"KO: {gene}")
    p16_indices = [i for i in range(len(all_seqs)) if p16_use[i]]
    if not p16_indices: p16_indices = list(range(min(100, len(all_seqs))))
    ko_dists = []
    for i in p16_indices:
        ko_seq = [t for t in all_seqs[i] if t != ko_token]
        ko_emb = get_emb(model, ko_seq, DEV)
        ko_dists.append(np.linalg.norm(ko_emb - p16neg_cent))
    ko_dist = np.mean(ko_dists)
    shift = orig_dist_pos - ko_dist
    results.append(dict(gene=gene, orig_dist=round(orig_dist_pos,4), ko_dist=round(ko_dist,4),
                       shift=round(shift,4), n_cells=len(p16_indices)))
    print(f"  orig: {orig_dist_pos:.4f} → KO: {ko_dist:.4f} | shift: {shift:+.4f}")

df = pd.DataFrame(results).sort_values("shift", ascending=False)
df.to_csv("/kaggle/working/mousegf_isp_results.csv", index=False)
print("\n=== MOUSE-GF ISP RESULTS ===")
print(df.to_string(index=False))
print(f"\nDONE ({time.time()-t0:.0f}s)")
