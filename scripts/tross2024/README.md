# Vendored: Tross et al. 2024 autoencoder code

Copied unmodified from https://github.com/mtross2/autoencoder_hyperspec_ref
at commit `24956e2908a487210ab23c7215dde61dd137edaa` on 2026-10-12. MIT
licence (`LICENSE`).

Paper: Tross et al. 2024, Plant Phenome J. 7:e20106, doi:10.1002/ppj2.20106.
Trained weights and full data: doi:10.6084/m9.figshare.24808491.v4 (not copied).

| File | Upstream path |
|---|---|
| `autoencoder_trainer/` | `autoencoder_trainer/` |
| `run_autoencoder_training.py` | `scripts/run_autoencoder_training.py` |
| `heritability_BroadSense.R` | `heritability_BroadSense.R` |
| `ae_architecture.json` | `trained_model/ae_architecture.json` |

Not copied: `requirements.txt` (conda export with local build paths, not
installable), `setup.py` and `README.md` (both name a nonexistent
`train_autoencoder.py`), `tests/`, `sample_data/`.

## Where code and paper differ

| | Paper | Code default | HypeSpec uses |
|---|---|---|---|
| Early-stop patience | 100 | 50 | 100 (`--patience 100`) |
| Train fraction | 5:1 | 0.8 | 5/6 (`--prop_split 0.8333`) |
| Dropout 0.3 after encoder layers 1–2 | not mentioned | yes | as code |

Input is raw plot-level reflectance, columns `X350`..`X2500`, a `PlotID`
column and a `trt` column, gzip CSV. No scaling. Run the training script with
this folder as the working directory (it adds the cwd to `sys.path`).

`heritability_BroadSense.R` needs replicated genotypes; the ZEAL BZ grid has
one plot per genotype, so it is kept for reference only.
