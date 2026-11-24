# Species ranges in China and worldwide have higher habitat quality yet face disproportionate decline

Data and code accompanying the manuscript submitted to *npj Biodiversity*. This package reproduces the statistical summaries, anthropogenic contribution indices, and main figures (Figs. 1–4) from archived 2000–2020 model outputs covering 1,820 IUCN Red List species across 96,649 nominal 10-km analysis grids in China and 29,671 terrestrial species (including 9,821 threatened species) worldwide.

## Package layout

- `china/` — China analysis scripts, InVEST parameters, 10-km model input table (`huigui_V2.parquet`), archived 2000–2020 SHAP attributions, grid keys and annual train/test splits, statistical test inputs/scripts, and source data for Figs. 1–3 and Supplementary Figs. 1–13.
- `global/` — Global extension scripts, cleaned configuration and land-cover crosswalks, species- and realm-level summary tables, and display grids/basemaps required to reproduce Fig. 4.
- `metadata/` — Variable definitions (`field_dictionary.csv`) and mapping from earlier internal filenames/table labels to the current manuscript (`legacy_name_crosswalk.csv`).
- `source_manifest.csv` — Complete file inventory with byte sizes and links to manuscript items.

*(Note: Code, configuration files, and metadata are mirrored on GitHub; large derived data tables and plotting grids are archived in this Figshare package.)*

## Quick reproduction

Run all commands from the root directory containing this `README.md` (keeping `china/` and `global/` together). Outputs are written to a sibling directory named `review-output` by default (`--output-dir` for Python scripts; first positional argument for R scripts).

**Environment.** Python 3.10–3.11 with `numpy`, `pandas`, `scipy` (>=1.11), `matplotlib`, `pyarrow`, and GDAL Python bindings (for example: `conda create -n hq-review -c conda-forge python=3.11 numpy pandas scipy matplotlib pyarrow gdal`). China figure scripts require R with `ggplot2`, `dplyr`, `grid`, `sf`, `ggmapcn`, `ragg`, `svglite`, `jsonlite`, and `boot`.

```shell
python -B china/scripts/summarize_shap.py
python -B china/scripts/anthropogenic_indices.py
python -B china/statistics/reproduce_statistics.py
python -B china/statistics/bootstrap_pairwise.py
Rscript china/statistics/bootstrap_table9.R
Rscript china/scripts/plot_figure1.R
Rscript china/scripts/plot_figure2.R
Rscript china/scripts/plot_figure3.R
python -B global/scripts/summarize_global_results.py
python -B global/scripts/plot_global_hq_signed_slopes.py
```

## China analysis

The manuscript attributions are derived from the January 20, 2026 archived SHAP matrix (`china/shap/`), which joins to grid keys (`grid_id`/`TID`) and projected coordinates (`china/model_inputs/`). Running `summarize_shap.py` reproduces all 252 rows (12 drivers × 21 years) of Supplementary Fig. 8, and `anthropogenic_indices.py` computes 2,029,629 grid-year values of ANCI and APCI from the negative and positive SHAP contributions of HFI, POP, GDP, and CE (normalized by the sum of absolute SHAP contributions across all 12 drivers plus `1e-6`), reproducing Fig. 3 and provincial means in Supplementary Table 14.

The 10-km model input table (`china/model_inputs/huigui_V2.parquet`) preserves all 236 numeric modeling fields and row order from the original workbook, paired with annual train/test splits (67,654 training and 28,995 test grids per year). The GWRF training and merging scripts document the January 29–30 model configuration; re-training local forests from scratch may produce slight stochastic differences relative to the archived January 20 attribution matrix. Supplementary Fig. 5 reports residual diagnostics from the separate annual prediction run.

**Grid-level trend archive (Fig. 1b).** `china/figure_data/Fig_1b_grid_trends.csv` is byte-identical to `Fig_3b_grid_trends.csv` referenced in the Supplementary Information (page S25) and the `Supplementary Data 1.xlsx` README. It provides fitted 2000–2020 annual habitat quality trends and uncertainty across all 96,649 10-km grids (p values are blank for 1,292 constant series).

Statistical reproduction scripts cover Supplementary Tables 4, 9, and 13 alongside the pairwise Welch (48 comparisons) and Mann–Whitney U (96 comparisons) tests in `Supplementary Data 1.xlsx`. The two bootstrap scripts reproduce the reported confidence intervals using 10,000 resamples, with the Kruskal–Wallis effect size standardized as `rank_eta_squared`.

## Global extension

Running `summarize_global_results.py` and `plot_global_hq_signed_slopes.py` reproduces the global statistics and Fig. 4 directly from `global/tables/` and `global/figure_data/`. Across 2000–2020, mean modeled habitat quality declined at 2.80× the fitted global-land rate across all 29,671 species and 3.72× across the 9,821 threatened species, while 5th-percentile habitat quality declined 182.1% and 133.4% faster than mean habitat quality, respectively. Across all species, 44.4% experienced declining 5th-percentile habitat quality, and 80.3% of those had a lower 5th-percentile slope than mean habitat quality slope. Low-quality habitat is defined separately by domain as HQ ≤ 0.239 in China and HQ < 0.2814 globally.

## Upstream data and reuse notes

All review commands above execute solely from the files distributed in this package. Re-running upstream spatial InVEST modeling from raw geospatial layers requires independently obtaining restricted or large third-party inputs listed in manuscript Table 1 (including IUCN Red List range polygons, annual land-cover rasters, and road/transport layers; see `global/config/parameters.md` for upstream configuration and China GIS entry points). Original IUCN species polygons and full-resolution source rasters are excluded in compliance with source data terms of use.
