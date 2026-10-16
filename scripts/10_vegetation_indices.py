#!/usr/bin/env python3
"""
Compute the 28 vegetation indices of Table 1 in "Phenotyping maize stay green
traits via in situ leaf hyperspectral reflectance sensing"
(https://www.sciencedirect.com/science/article/pii/S2772375526001875) for every
BZ plot, from the splice-corrected plot-mean reflectance (1-nm bands).

Equations are as printed in that table, except two evident typos (noted in the
definitions sheet): NDRE's denominator (printed R780 + R730) uses R790 + R720,
and Chl NDI's denominator (printed R750 + 705) uses R750 + R705.

In : output/bz_plot_reflectance.csv.gz, output/bz_plot_meta.csv,
     output/outlier_group_plots.csv
Out: output/vegetation_indices_table1.xlsx  (sheets: indices, definitions)
     output/vegetation_indices_table1.csv   (indices sheet only)
Run: .venv/Scripts/python scripts/10_vegetation_indices.py
"""
import pandas as pd

spec = pd.read_csv("output/bz_plot_reflectance.csv.gz")
meta = pd.read_csv("output/bz_plot_meta.csv")
sen = set(pd.read_csv("output/outlier_group_plots.csv").PlotID)


def R(nm):
    return spec[f"X{nm}"]


swir = spec[[f"X{w}" for w in range(1500, 1751)]]
# (group, acronym, name, equation as used, reference in Table 1, note)
defs = [
    ("Chlorophyll content and greenness", "CIgreen", "Green Chlorophyll Index", "R730/R530 - 1", "[35]", ""),
    ("Chlorophyll content and greenness", "RARSa", "Ratio Analysis of Reflectance Spectra Chlorophyll a", "R675/R700", "[42]", ""),
    ("Chlorophyll content and greenness", "RARSb", "Ratio Analysis of Reflectance Spectra Chlorophyll b", "R675/(R650 x R700)", "[43]", ""),
    ("Chlorophyll content and greenness", "NDVI", "Normalized Difference Vegetation Index", "(R800 - R680)/(R800 + R680)", "[44]", ""),
    ("Chlorophyll content and greenness", "OSAVI", "Optimized Soil Adjusted Vegetation Index", "(1 + 0.16)(R800 - R670)/(R800 + R670 + 0.16)", "[45]", ""),
    ("Chlorophyll content and greenness", "NDRE", "Normalized Difference Red Edge Index", "(R790 - R720)/(R790 + R720)", "[46]", "Table prints denominator R780 + R730; standard symmetric form used"),
    ("Chlorophyll content and greenness", "PRI", "Photochemical reflectance index", "(R531 - R570)/(R531 + R570)", "[47]", ""),
    ("Chlorophyll content and greenness", "REIP", "Red-edge inflection point", "700 + 40 x [((R670 + R780)/2 - R700)/(R740 - R700)]", "[48]", "Units: nm"),
    ("Chlorophyll content and greenness", "NNDVI", "Narrow-band Normalized Difference Vegetation Index", "(R775 - R670)/(R775 + R670)", "[49]", ""),
    ("Chlorophyll content and greenness", "CI", "Chlorophyll index", "R880/R590 - 1", "[50]", ""),
    ("Chlorophyll content and greenness", "CIred-edge", "Chlorophyll index (red edge)", "R850/R730 - 1", "[50]", ""),
    ("Chlorophyll content and greenness", "DD", "Double difference index", "(R749 - R720) - (R701 - R672)", "[51]", ""),
    ("Chlorophyll content and greenness", "SIPI", "Structure intensive pigment index", "(R800 - R445)/(R800 - R680)", "[52]", ""),
    ("Chlorophyll content and greenness", "Chl NDI", "Normalized Difference 750/705", "(R750 - R705)/(R750 + R705)", "[53]", "Table prints denominator R750 + 705; R705 used"),
    ("Chlorophyll content and greenness", "Red-edge", "Simple ratio", "R750/R700", "[54]", ""),
    ("Chlorophyll content and greenness", "GRE", "Greenness", "500 x (R800/R710 - 1) + 27", "[35]", ""),
    ("Chlorophyll content and greenness", "NDVI2", "Normalized Difference Vegetation Index (760/710)", "(R760 - R710)/(R760 + R710)", "source study", ""),
    ("Senescence", "PSRI", "Plant Senescence Reflectance Index", "(R678 - R500)/R750", "[55]", ""),
    ("Senescence", "NDIS_2b", "Normalized Difference Index for Senescence (2 bands)", "(R2190 - R1510)/(R2190 + R1510)", "source study", ""),
    ("Senescence", "RDIS_3b", "Ratio of Difference Index for Senescence (3 bands)", "(R1280 - R1250)/(R1250 - R1180)", "source study", ""),
    ("Water status", "NDWI_1", "Normalized Difference Water Index", "(R864 - R1610)/(R864 + R1610)", "[56]", ""),
    ("Water status", "NDWI_2", "Normalized Difference Water Index", "(R860 - R1241)/(R860 + R1241)", "[57]", ""),
    ("Water status", "NDWI_3", "Normalized Difference Water Index", "(R864 - R2200)/(R864 + R2200)", "[56]", ""),
    ("Water status", "NDWI_4", "Normalized Difference Water Index", "(R559 - R864)/(R559 + R864)", "[58]", ""),
    ("Water status", "WI", "Water Index", "R900/R970", "[59]", ""),
    ("Water status", "WI/NDVI", "Ratio WI/NDVI", "(R900/R970)/[(R800 - R680)/(R800 + R680)]", "[59]", ""),
    ("Water status", "MDWI", "Maximum difference water index", "(Rmax - Rmin)/(Rmax + Rmin), R over 1500-1750 nm", "[60]", ""),
    ("Water status", "MSI", "Moisture Stress Index", "R1600/R820", "[61]", ""),
]

ndvi = (R(800) - R(680)) / (R(800) + R(680))
v = {
    "CIgreen": R(730) / R(530) - 1,
    "RARSa": R(675) / R(700),
    "RARSb": R(675) / (R(650) * R(700)),
    "NDVI": ndvi,
    "OSAVI": 1.16 * (R(800) - R(670)) / (R(800) + R(670) + 0.16),
    "NDRE": (R(790) - R(720)) / (R(790) + R(720)),
    "PRI": (R(531) - R(570)) / (R(531) + R(570)),
    "REIP": 700 + 40 * (((R(670) + R(780)) / 2 - R(700)) / (R(740) - R(700))),
    "NNDVI": (R(775) - R(670)) / (R(775) + R(670)),
    "CI": R(880) / R(590) - 1,
    "CIred-edge": R(850) / R(730) - 1,
    "DD": (R(749) - R(720)) - (R(701) - R(672)),
    "SIPI": (R(800) - R(445)) / (R(800) - R(680)),
    "Chl NDI": (R(750) - R(705)) / (R(750) + R(705)),
    "Red-edge": R(750) / R(700),
    "GRE": 500 * (R(800) / R(710) - 1) + 27,
    "NDVI2": (R(760) - R(710)) / (R(760) + R(710)),
    "PSRI": (R(678) - R(500)) / R(750),
    "NDIS_2b": (R(2190) - R(1510)) / (R(2190) + R(1510)),
    "RDIS_3b": (R(1280) - R(1250)) / (R(1250) - R(1180)),
    "NDWI_1": (R(864) - R(1610)) / (R(864) + R(1610)),
    "NDWI_2": (R(860) - R(1241)) / (R(860) + R(1241)),
    "NDWI_3": (R(864) - R(2200)) / (R(864) + R(2200)),
    "NDWI_4": (R(559) - R(864)) / (R(559) + R(864)),
    "WI": R(900) / R(970),
    "WI/NDVI": (R(900) / R(970)) / ndvi,
    "MDWI": (swir.max(axis=1) - swir.min(axis=1)) / (swir.max(axis=1) + swir.min(axis=1)),
    "MSI": R(1600) / R(820),
}
assert [d[1] for d in defs] == list(v)

m = meta.set_index("PlotID").loc[spec.PlotID]
out = pd.DataFrame({
    "PlotID": spec.PlotID.values,
    "scan_date": m["date"].values,
    "Row": m["Row"].values,
    "Range": m["Range"].values,
    "Genotype": m["Genotype"].values,
    "Taxa": m["Taxa"].values,
    "taxa_mex": m["taxa_mex"].values,
    "senesced": spec.PlotID.isin(sen).values,
    "Predicted_N_mean": m["predN_mean"].values,
})
for k, s in v.items():
    out[k] = s.values

dd = pd.DataFrame(defs, columns=["group", "acronym", "name", "equation_used", "reference_in_table1", "note"])
dd.insert(4, "source", "Table 1, 'Phenotyping maize stay green traits via in situ leaf hyperspectral "
          "reflectance sensing', https://www.sciencedirect.com/science/article/pii/S2772375526001875")

with pd.ExcelWriter("output/vegetation_indices_table1.xlsx") as xw:
    out.to_excel(xw, sheet_name="indices", index=False)
    dd.to_excel(xw, sheet_name="definitions", index=False)
out.to_csv("output/vegetation_indices_table1.csv", index=False)

g = out[~out.senesced]
summ = pd.DataFrame({
    "median_all": out[list(v)].median(),
    "median_nonsenesced": g[list(v)].median(),
    "median_senesced": out.loc[out.senesced, list(v)].median(),
    "rho_predN_nonsenesced": g[list(v)].corrwith(g.Predicted_N_mean, method="spearman"),
    "n_nonfinite": out[list(v)].apply(lambda c: (~pd.Series(c).replace([float("inf"), -float("inf")], float("nan")).notna()).sum()),
}).round(3)
print(summ.to_string())
