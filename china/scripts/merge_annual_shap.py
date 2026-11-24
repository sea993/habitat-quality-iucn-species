"""Join all annual refit exports by projected coordinate to the study grid."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from shap_data import CHINA, DRIVERS, YEARS, column


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = CHINA.parent.parent / 'review-output/refit'
    parser.add_argument('--input-dir', type=Path, default=root / 'annual')
    parser.add_argument('--output-dir', type=Path, default=root)
    args = parser.parse_args()
    keys = pd.read_csv(CHINA/'model_inputs/grid_key.csv').set_index('grid_id',verify_integrity=True)
    tree = cKDTree(keys[['X','Y']])
    frames = []
    for year in YEARS:
        data = pd.read_csv(args.input_dir / f'GWRF_SHAP_Result_{year}.csv')
        distance, ix = tree.query(data[['X','Y']])
        if len(ix)!=len(keys) or len(np.unique(ix))!=len(keys) or distance.max()>.01:
            raise ValueError(f'{year}: coordinate match failed (tolerance 0.01 m)')
        data.index = keys.index[ix]
        names={'SHAP_'+(d if d in ('DEM','SLOPE') else d+str(year)):column(d,year) for d in DRIVERS}
        frame=data.loc[keys.index,list(names)].rename(columns=names)
        if not np.isfinite(frame.to_numpy()).all():
            raise ValueError(f'{year}: incomplete SHAP values')
        frames.append(frame)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    pd.concat(frames,axis=1).to_csv(args.output_dir/'refit_grid_SHAP_2000_2020.csv')


if __name__ == '__main__':
    main()
