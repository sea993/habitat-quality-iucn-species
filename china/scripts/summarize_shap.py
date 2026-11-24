"""Recompute the 21 by 12 driver summary for Supplementary Fig. 8."""
import argparse
from pathlib import Path
import pandas as pd
from shap_data import DRIVERS, YEARS, OUTPUT, column, load_grid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    out = parser.parse_args().output_dir
    _, grid = load_grid()
    rows = []
    for year in YEARS:
        columns = [column(driver, year) for driver in DRIVERS]
        signed = grid[columns].mean()
        magnitudes = grid[columns].abs().mean()
        for driver, name in zip(DRIVERS, columns):
            rows.append(dict(Year=year, Driver=driver, SHAP_column='SHAP_' + name,
                             Mean_signed_SHAP=signed[name], Mean_absolute_SHAP=magnitudes[name],
                             Annual_absolute_contribution_percent=100*magnitudes[name]/magnitudes.sum(),
                             N_grid_cells=len(grid)))
    summary = pd.DataFrame(rows)
    ranks = summary.groupby('Driver').Mean_absolute_SHAP.mean().rank(ascending=False).astype(int)
    summary['Rank_by_2000_2020_mean_absolute_SHAP'] = summary.Driver.map(ranks)
    out.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out / 'Figure_S8_source.csv', index=False)


if __name__ == '__main__':
    main()
