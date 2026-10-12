# HypeSpec
For analyzing hyperspectral data taken on maize leaves.

Unsupervised leaf-reflectance phenotypes for the ZEAL (BZea) teosinte
introgression population, following Tross et al. 2024 (Plant Phenome J.
7:e20106, doi:10.1002/ppj2.20106). Data: FieldSpec4 scans of leaf 2 above the
ear, CLY25-B5 block 2 (921 plots, 5 scans per plot, 2025-07-18 to 07-26).

## Pipeline

```mermaid
flowchart TD
    A["Raw scans (.asd)<br/>4,615 FieldSpec4 scans, 350–2500 nm"] --> B["Indico Pro: Asd to Spc<br/>log(1/R) per scan (.spc)"]
    B --> C["01_plot_spectra.R<br/>R = 10^−log(1/R) · splice correction at 1000 / 1800 nm<br/>mean of 5 scans per plot · reflectance 0–1 check"]
    M["Collection workbook (BZ-remap)<br/>predicted-N plot key"] --> C
    C --> D[("bz_plot_reflectance.csv.gz<br/>921 plots × 2151 bands")]
    C --> E[("bz_plot_meta.csv<br/>genotype, taxa, notes,<br/>Predicted_N, M-distance")]

    D --> P["02_pca.R<br/>PCA, 10 PCs"]
    D --> T["07_tross_pretrained.py<br/>Tross trained encoder, 10 LVs"]
    D --> R["06_train_autoencoder.py<br/>autoencoder retrained on ZEAL, 10 LVs"]

    P --> P2["03–05: PCA figures and<br/>component table"]
    E --> P2

    P --> H["SpATS BLUEs and H² per component<br/>(planned)"]
    T --> H
    R --> H
    E --> H
    H --> Q["R/qtl scans and peaks in ZeaL_QTL<br/>(planned)"]
    E --> G["Annotation: Spearman ρ with Predicted_N<br/>(planned for LVs)"]
    T --> G
    R --> G

    classDef planned stroke-dasharray: 5 5
    class H,Q,G planned
```

## Layout

| Folder | Contents |
|---|---|
| `scripts/` | Numbered pipeline scripts; `tross2024/` holds the vendored Tross code (see its README) |
| `data/` | Inputs (gitignored) |
| `output/` | Results (gitignored) |

R scripts run from the repo root with `Rscript scripts/<name>.R`. Python
scripts use the local environment: `.venv/Scripts/python scripts/<name>.py`
(TensorFlow, tf_keras, pandas, scikit-learn, h5py).
