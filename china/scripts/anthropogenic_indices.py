"""Compute the two signed anthropogenic contribution indices for every grid-year."""
import argparse
from pathlib import Path
from shap_data import DRIVERS, YEARS, OUTPUT, column, load_grid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    out = parser.parse_args().output_dir
    out.mkdir(parents=True, exist_ok=True)
    keys, grid = load_grid()
    destination = out / 'grid_year_indices.csv'
    temporary = destination.with_suffix('.partial')
    for year in YEARS:
        total = grid[[column(d, year) for d in DRIVERS]].abs().sum(axis=1)
        human = grid[[column(d, year) for d in ('HFI', 'POP', 'GDP', 'CE')]]
        # The archived calculation adds 1e-6 for every cell to prevent division by zero.
        rows = keys.assign(year=year, ANCI=-100*human.clip(upper=0).sum(axis=1)/(total+1e-6),
                           APCI=100*human.clip(lower=0).sum(axis=1)/(total+1e-6))
        rows.to_csv(temporary, mode='w' if year == 2000 else 'a', header=year == 2000)
    temporary.replace(destination)


if __name__ == '__main__':
    main()
