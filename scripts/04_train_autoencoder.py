#!/usr/bin/env python3
"""
Train the Tross et al. 2024 autoencoder on the BZ plot-level spectra.

Runs the vendored, unmodified trainer (scripts/tross2024) with the paper's
settings where paper and code disagree (patience 100, 5:1 train:validation).
Everything else is the vendored code's default: 1000 epochs max, SGD lr 0.1,
MAE loss, split seed 0.

Two shims, neither changes the method:
  * the vendored code was written for Keras 2 (TF 2.4-2.8); modern TF ships
    Keras 3, so Keras 2 (tf_keras) is used and aliased as `keras`.
  * TF weight initialisation is seeded (Tross seeded only the split) so a
    rerun gives the same model.

In : output/bz_plot_reflectance.csv.gz
Out: output/ae_<N>lv/  (LVs_Reflectance_Data_<tag>.csv, losses, weights, model)
Run: .venv/Scripts/python scripts/04_train_autoencoder.py [--num_variables 10]
"""
import argparse
import os
import sys

os.environ["TF_USE_LEGACY_KERAS"] = "1"  # tensorflow.keras -> tf_keras (Keras 2)
import tf_keras  # noqa: E402

sys.modules["keras"] = tf_keras  # vendored code does `from keras.callbacks import ...`
sys.modules["keras.callbacks"] = tf_keras.callbacks
import tensorflow as tf  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "tross2024"))
from autoencoder_trainer import autoencoder_trainer  # noqa: E402

p = argparse.ArgumentParser()
p.add_argument("--dataset", default="output/bz_plot_reflectance.csv.gz")
p.add_argument("--num_variables", type=int, default=10)
p.add_argument("--epochs", type=int, default=1000)
p.add_argument("--patience", type=int, default=100)  # paper; code default 50
p.add_argument("--prop_split", type=float, default=5 / 6)  # paper 5:1; code default 0.8
p.add_argument("--learning_rate", type=float, default=0.1)
p.add_argument("--random_seed", type=int, default=0)
p.add_argument("--exclude", help="CSV with a PlotID column; those plots are dropped")
p.add_argument("--tag", default="", help="suffix for the output folder and file tag, e.g. _green")
a = p.parse_args()

tf_keras.utils.set_random_seed(a.random_seed)
out = f"output/ae_{a.num_variables}lv{a.tag}"
os.makedirs(out, exist_ok=True)
dataset = a.dataset
if a.exclude:
    import pandas as pd

    d = pd.read_csv(a.dataset)
    d = d[~d.PlotID.isin(pd.read_csv(a.exclude).PlotID)]
    dataset = os.path.join(out, "input_spectra.csv.gz")
    d.to_csv(dataset, index=False, compression="gzip")
    print(f"exclude: {len(d)} plots remain -> {dataset}", flush=True)
trainer = autoencoder_trainer.AutoencoderTrainer(
    dataset_path=dataset, codings_size=a.num_variables, epochs=a.epochs,
    patience=a.patience, learning_rate=a.learning_rate, random_seed=a.random_seed,
    prop_split=a.prop_split, output_loc=out, ts=f"bz_{a.num_variables}lv{a.tag}",
)
print(f"train {len(trainer.df_train)} / validation {len(trainer.df_validation)} plots", flush=True)
trainer.train()
