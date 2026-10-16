#!/usr/bin/env python3
"""
Apply the Tross et al. 2024 trained encoder to the BZ plot-level spectra.

The encoder (2151 -> 2200 -> 3000 -> 2024 -> 10; SELU, linear latent layer;
dropout is inactive at inference) is evaluated directly from the published
weights with numpy, so no Keras version issues arise.

1. VERIFY: encode Tross's own published spectra and compare with their
   published latent variables. Rows are matched within GenotypeID (two reps,
   and the files differ in row count), taking the closest computed row.
2. APPLY: encode our spectra (same 350-2500 nm, 1 nm, raw reflectance input).

In : data/tross2024/weights.h5                       (figshare 24808491 v4)
     data/tross2024/maize_2020_unl_reflectance_traits.csv
     data/tross2024/maize_2020_unl_latentVariables.csv
     output/bz_plot_reflectance.csv.gz
Out: output/tross_pretrained_lvs.csv   (PlotID, LV1..LV10)
Run: .venv/Scripts/python scripts/05_tross_pretrained.py
"""
import h5py
import numpy as np
import pandas as pd

W = "data/tross2024/weights.h5"
ALPHA, SCALE = 1.6732632423543772, 1.0507009873554805  # SELU constants


def selu(x):
    return SCALE * np.where(x > 0, x, ALPHA * np.expm1(np.minimum(x, 0)))


with h5py.File(W, "r") as f:
    enc = f["model_weights"][[k for k in f["model_weights"] if len(f["model_weights"][k]) == 4][0]]
    layers = sorted(enc, key=lambda n: int(n.split("_")[1]))  # dense_291..294
    params = [(enc[l]["kernel:0"][()], enc[l]["bias:0"][()]) for l in layers]
assert [k.shape for k, _ in params] == [(2151, 2200), (2200, 3000), (3000, 2024), (2024, 10)]


def encode(X):
    h = X.astype(np.float32)
    for i, (k, b) in enumerate(params):
        h = h @ k + b
        if i < len(params) - 1:
            h = selu(h)
    return h


lvcols = [f"LV{i}" for i in range(1, 11)]
xcols = [f"X{w}" for w in range(350, 2501)]

# --- 1. verify against published latent variables ---------------------------
ref = pd.read_csv("data/tross2024/maize_2020_unl_reflectance_traits.csv")
pub = pd.read_csv("data/tross2024/maize_2020_unl_latentVariables.csv")
ok = ref[xcols].notna().all(axis=1)
mine = pd.DataFrame(encode(ref.loc[ok, xcols].to_numpy()), columns=lvcols)
mine["GenotypeID"] = ref.loc[ok, "GenotypeID"].to_numpy()
diffs = []
for g, rows in pub.groupby("GenotypeID"):
    cand = mine.loc[mine.GenotypeID == g, lvcols].to_numpy()
    if len(cand) == 0:
        continue
    for r in rows[lvcols].to_numpy():
        diffs.append(np.abs(cand - r).max(axis=1).min())
diffs = np.array(diffs)
print(f"verify: {len(diffs)} of {len(pub)} published rows matched; "
      f"max |diff| median {np.median(diffs):.2e}, 99th pct {np.quantile(diffs, 0.99):.2e}, "
      f"max {diffs.max():.2e}; LV range in published data "
      f"{pub[lvcols].to_numpy().min():.2f} to {pub[lvcols].to_numpy().max():.2f}")

# --- 2. apply to BZ ----------------------------------------------------------
bz = pd.read_csv("output/bz_plot_reflectance.csv.gz")
out = pd.DataFrame(encode(bz[xcols].to_numpy()), columns=lvcols)
out.insert(0, "PlotID", bz.PlotID)
out.to_csv("output/tross_pretrained_lvs.csv", index=False)
print(f"wrote output/tross_pretrained_lvs.csv ({len(out)} plots)")
print(out[lvcols].describe().loc[["mean", "std", "min", "max"]].round(3).to_string())
