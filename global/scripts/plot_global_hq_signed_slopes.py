"""Global HQ changes, with a shared baseline year and a continuous ochre–cyan map."""

import json

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
from matplotlib.transforms import Bbox
from osgeo import gdal, ogr

from pathlib import Path
import argparse
from figure4_layout import (COLORS, GROUPS, HEIGHT, WIDTH, axis_mm, map_panel, new_figure,
                            paired_panels, plt, small_title, _fit_gamut, _hex_to_rgb,
                            _linear_to_srgb, _oklab_to_linear, _rgb_to_oklab)

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / 'tables'
STEM = 'Fig_4'
# Axis transition in plotted units of 0.001 HQ/year.
LINEAR_WIDTH = 0.5
YEARS = np.arange(2000, 2021)
NODES = (('#F2AA87', '#D67852', '#B75039'),
         ('#FAEBCB', '#ECE5C7', '#C0D5BE'),
         ('#A2E3E0', '#5ABEC8', '#258DAA'))
NO_DATA = '#EEEEEF'


def habitat_colors():
    """Interpolate the nine design anchors in OKLab on two continuous axes."""
    coords = np.linspace(0, 1, 257)
    knots = [0, .5, 1]
    nodes = _rgb_to_oklab(np.array([[_hex_to_rgb(c) for c in row] for row in NODES]))
    rows = np.array([np.column_stack([np.interp(coords, knots, channel) for channel in row.T])
                     for row in nodes])
    magnitude = np.abs(2 * coords - 1)[:, None, None]
    # Gentle color progression near zero keeps modest changes legible at print size.
    strength = magnitude * (1.5 - .5 * magnitude)
    ends = np.where((coords < .5)[:, None, None], rows[0], rows[2])
    lab = rows[1] + strength * (ends - rows[1])
    return np.rint(255 * _linear_to_srgb(_oklab_to_linear(_fit_gamut(lab)))).astype(np.uint8)


def land_data():
    annual = pd.read_csv(TABLE / 'background_annual.csv')
    cfg = json.loads((TABLE / 'global_low_quality_threshold.json').read_text(encoding='utf-8'))
    frames = {}
    with np.load(TABLE / 'background_histograms.npz') as source:
        years = source['years']
        if not np.array_equal(years, YEARS):
            raise ValueError('Missing land histogram years')
        for scope in ('Global',):
            table = annual[annual.scope == scope].set_index('year').sort_index()
            if not np.array_equal(table.index, YEARS):
                raise ValueError(f'Unexpected years for {scope}')
            hist = source['hist'][list(source['scopes']).index(scope)]
            counts = hist.sum(axis=1)
            if not np.array_equal(counts, table.grid_cells) or table.grid_cells.nunique() != 1:
                raise ValueError(f'Common land area mismatch for {scope}')
            idx = (np.cumsum(hist, axis=1) >= counts[:, None] * .05).argmax(axis=1)
            p5 = (idx + 1) * cfg['histogram_bin_width']
            if not np.allclose(p5, table.q05_hq_histogram_upper_edge, rtol=0, atol=1e-12):
                raise ValueError(f'P5 mismatch for {scope}')
            frames[scope] = table
    return frames


def change(series):
    values = series.to_numpy()
    return values - values[0]


def trajectory_row(fig, species, land):
    specs = [('b', 'mean_hq', 'mean_hq', 'HQ change since 2000', .002, '%.3f'),
             ('c', 'q05_hq', None, 'HQ change since 2000', .005, '%.3f'),
             ('d', 'low_quality_pct', 'low_quality_pct', 'Change (percentage points)', .25, '%.2f')]
    axes = []
    for i, (letter, metric, field, ylabel, tick, fmt) in enumerate(specs):
        ax = axis_mm(fig, (15 + i * 60, 58, 42, 28))
        for group in GROUPS:
            data = species[(species.metric == metric) & (species.group == group)].set_index('year').sort_index()
            if not np.array_equal(data.index, YEARS):
                raise ValueError(f'Unexpected trajectory years for {group} {metric}')
            if not np.allclose(change(data.estimate), data.change_from_2000):
                raise ValueError(f'Unexpected reference year for {group} {metric}')
            ax.fill_between(YEARS, data.change_ci95_low.to_numpy(), data.change_ci95_high.to_numpy(),
                            color=COLORS[group], alpha=.17, lw=0)
            ax.plot(YEARS, data.change_from_2000, color=COLORS[group], lw=1.15)
        if field:
            global_change = change(land['Global'][field])
            ax.plot(YEARS, global_change, color=COLORS['Global land'], ls=(0, (4, 2)), lw=.9)
        ax.set(xlim=(1999.25, 2020.75), xticks=[2000, 2010, 2020], ylabel=ylabel)
        ax.margins(y=.09)
        ax.yaxis.set_major_locator(MultipleLocator(tick))
        ax.yaxis.set_major_formatter(FormatStrFormatter(fmt))
        small_title(ax, letter, '')
        axes.append(ax)
    handles = [Line2D([], [], color=COLORS[g], lw=1.15, label=g) for g in GROUPS]
    handles.append(Line2D([], [], color=COLORS['Global land'], lw=.9, ls=(0, (4, 2)), label='Global land'))
    fig.legend(handles=handles, loc='center left', ncol=3, borderaxespad=0, handlelength=2.0,
               columnspacing=1.4, bbox_to_anchor=(15/WIDTH, 94/HEIGHT))
    return axes


def comparison_row(fig, slopes, effects):
    pair, realms = paired_panels(fig, slopes, effects, row=(20, 25.5), bar_box=(27, 8.5, 30, 1.4),
                                 panel_title=small_title)
    old_mesh = pair.collections[0]
    grey = old_mesh.get_cmap()
    old_mesh.colorbar.remove()
    old_mesh.remove()
    pair.lines[0].remove()
    xy = slopes[['mean_hq_slope','q05_hq_slope']].to_numpy() * 1000
    transformed = np.arcsinh(xy / LINEAR_WIDTH)
    limit = np.abs(transformed).max() * 1.025
    mesh = pair.hexbin(*transformed.T, extent=(-limit, limit, -limit, limit),
                       gridsize=(46, 25), mincnt=1, bins='log', cmap=grey,
                       linewidths=0, rasterized=True)
    if int(mesh.get_array().sum()) != len(slopes):
        raise ValueError('Paired slope display did not retain every species')
    pair.plot([-limit, limit], [-limit, limit], ls='--', color='#666666', lw=.65, zorder=0)
    ticks = np.array([-30,-1,0,1,30])
    for axis in (pair.xaxis, pair.yaxis):
        axis.set_ticks(np.arcsinh(ticks / LINEAR_WIDTH),
                       labels=[str(int(v)).replace('-', '−') for v in ticks])
    pair.set(xlim=(-limit,limit), ylim=(-limit,limit))
    bar = fig.colorbar(mesh, cax=axis_mm(fig,(61,20,1.8,25.5)), orientation='vertical')
    count_ticks = [10**i for i in range(1 + int(np.floor(np.log10(mesh.get_array().max()))))]
    bar.set_ticks(count_ticks, labels=[f'{n:,}' for n in count_ticks])
    bar.ax.minorticks_off()
    bar.ax.tick_params(length=1.6,pad=1.4)
    bar.set_label('Species per hexagon', labelpad=2)
    bar.outline.set_visible(False)
    pair.set_xlabel('Mean HQ slope (0.001/yr)')
    pair.set_ylabel('5th-percentile HQ slope\n(0.001/yr)')
    # Widen the gutter between the right-hand axis of panel e and the realm names.
    realms.set_position((98/WIDTH, 20/HEIGHT, 79/WIDTH, 25.5/HEIGHT))
    realms.set_xlabel('5th-percentile − mean HQ slope (0.001/yr)')
    realms.spines['left'].set_visible(False)
    realms.set_axisbelow(True)
    realms.grid(axis='y', lw=.32, color='#e9e9e9')
    return pair, realms


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT.parent.parent / 'review-output' / 'global')
    FIG = parser.parse_args().output_dir.resolve()
    FIG.mkdir(parents=True, exist_ok=True)
    gdal.UseExceptions()
    ogr.UseExceptions()
    species = pd.read_csv(TABLE / 'group_annual_trajectories.csv').query('scope == "Global"')
    effects = pd.read_csv(TABLE / 'group_effects_ci.csv')
    slopes = pd.read_csv(TABLE / 'species_slopes.csv').query('scope == "Global"')
    land = land_data()
    fig = new_figure()
    world, metadata = map_panel(fig, habitat_colors(), NO_DATA, baseline_table=True)
    starts = species.query("metric == 'mean_hq' and year == 2000").set_index('group').estimate
    left, top = metadata['baseline_anchor_mm']
    dots = axis_mm(fig, (left+1.2, top-8.2, 18.8, 6.6))
    dots.patch.set_visible(False)
    values = [starts[group] for group in GROUPS] + [land['Global'].mean_hq.iloc[0]]
    for row, (group, value) in enumerate(zip((*GROUPS, 'Global land'), values)):
        dots.plot(value, 2-row, 'o', ms=3.4, color=COLORS[group], zorder=3)
    dots.set(xlim=(.575,.825), ylim=(-.5,2.5), yticks=[], xticks=[.6,.7,.8])
    dots.spines['left'].set_visible(False)
    dots.spines['bottom'].set_linewidth(.4)
    dots.tick_params(axis='x', length=1.5, pad=1)
    dots.set_axisbelow(True)
    dots.grid(axis='x', color='#ddd8cf', linewidth=.4)
    fig.text(left/WIDTH, top/HEIGHT, 'Mean HQ in 2000', size=7, weight='bold', va='center')
    series = trajectory_row(fig, species, land)
    paired = comparison_row(fig, slopes, effects)
    # Use the actual painted pixels to remove empty canvas consistently across formats.
    fig.canvas.draw()
    painted = np.min(np.asarray(fig.canvas.buffer_rgba())[...,:3], axis=2) < 242
    iy, ix = np.nonzero(painted)
    height = painted.shape[0]
    trim = Bbox.from_extents(ix.min()/fig.dpi, (height-1-iy.max())/fig.dpi,
                             (ix.max()+1)/fig.dpi, (height-iy.min())/fig.dpi).padded(.04)
    fig.savefig(FIG / f'{STEM}.pdf', dpi=600, bbox_inches=trim)
    fig.savefig(FIG / f'{STEM}.svg', dpi=600, bbox_inches=trim)
    fig.savefig(FIG / f'{STEM}.png', dpi=300, bbox_inches=trim)
    plt.close(fig)
    metadata.update(color_nodes=NODES, no_data=NO_DATA, change_reference_year=2000)
    (FIG / f'{STEM}_display.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')



if __name__ == '__main__':
    main()
