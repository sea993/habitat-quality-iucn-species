"""Read the January 20 manuscript grid, checking identities before joining coordinates."""
from pathlib import Path
import numpy as np
import pandas as pd

CHINA = Path(__file__).resolve().parents[1]
OUTPUT = CHINA.parent.parent / 'review-output' / 'china'
YEARS = range(2000, 2021)
DRIVERS = 'NDVI GDP PRE POP CE HFI TMP CSDI GSL WSDI DEM SLOPE'.split()


def column(driver, year):
    return f'{driver}_{year}' if driver in ('DEM', 'SLOPE') else f'{driver}{year}'


def load_grid():
    keys = pd.read_csv(CHINA / 'model_inputs/grid_key.csv').set_index('grid_id', verify_integrity=True)
    values = pd.read_csv(CHINA / 'shap/Manuscript_grid_SHAP_2000_2020.csv').set_index('grid_id', verify_integrity=True)
    expected = [column(d, y) for y in YEARS for d in DRIVERS]
    if (len(keys) != 96649 or not keys.TID.is_unique or keys.isna().any().any()
            or set(values.columns) != set(expected) or set(values.index) != set(keys.index)):
        raise ValueError('Grid identities, coordinates or driver-year columns are incomplete')
    if not np.isfinite(values.to_numpy()).all():
        raise ValueError('Nonfinite SHAP value in manuscript grid')
    return keys, values.loc[keys.index]
