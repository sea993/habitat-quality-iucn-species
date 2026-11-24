"""Verify published global results from the annual species records and common-land baseline."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tables', type=Path, default=ROOT / 'tables')
    parser.add_argument('--output-dir', type=Path, default=ROOT.parent.parent / 'review-output/global')
    args = parser.parse_args()
    columns = ['scope', 'species_id', 'category', 'year', 'mean_hq', 'q05_hq']
    parts = [chunk.loc[chunk.scope == 'Global'] for chunk in pd.read_csv(
        args.tables / 'species_annual_metrics.csv', usecols=columns, chunksize=150000, float_precision='round_trip')]
    annual = pd.concat(parts)
    years = np.arange(2000, 2021)
    x = years - years.mean()
    slope = {}
    for metric in ('mean_hq', 'q05_hq'):
        values = annual.pivot(index='species_id', columns='year', values=metric).loc[:, years]
        assert values.shape == (29671, 21) and np.isfinite(values.to_numpy()).all()
        matrix = values.to_numpy()
        slope[metric] = (matrix - matrix[:, :1]) @ x / (x @ x)
    df = pd.DataFrame(slope, index=values.index)
    ref = pd.read_csv(args.tables / 'species_slopes.csv', float_precision='round_trip').query("scope == 'Global'").set_index('species_id').loc[df.index]
    np.testing.assert_allclose(df, ref[['mean_hq_slope', 'q05_hq_slope']], atol=1e-14, rtol=0)
    land = pd.read_csv(args.tables / 'background_annual.csv').query("scope == 'Global'").set_index('year').loc[years]
    baseline = np.dot(land.mean_hq - land.mean_hq.iloc[0], x) / (x @ x)
    groups = [('all', df), ('threatened', df.loc[ref.category.isin(['VU','EN','CR'])])]
    assert len(groups[1][1]) == 9821
    results = {}
    for name, group in groups:
        results[name + '_mean_to_land'] = group.mean_hq.mean() / baseline
        results[name + '_P5_excess_percent'] = 100*(group.q05_hq.mean()/group.mean_hq.mean()-1)
        declining = group.loc[group.q05_hq < 0]
        results[name + '_P5_declines_percent'] = 100*len(declining)/len(group)
        results[name + '_faster_among_decliners_percent'] = 100*(declining.q05_hq < declining.mean_hq).mean()
    expected = [2.80, 182.1, 44.4, 80.3, 3.72, 133.4, 38.6, 80.2]
    for (key, value), published in zip(results.items(), expected):
        digits = 2 if key.endswith('to_land') else 1
        if round(value, digits) != published:
            raise ValueError(f'{key}: {value} differs from manuscript {published}')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report = json.dumps(results, indent=2)
    (args.output_dir / 'global_results.json').write_text(report+'\n', encoding='utf8')
    print(report)


if __name__ == '__main__':
    main()
