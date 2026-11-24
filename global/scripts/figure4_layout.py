"""Map projection, colour interpolation and panel layout for Fig. 4."""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.ticker import FuncFormatter, FormatStrFormatter, MaxNLocator, MultipleLocator
from matplotlib.patches import Polygon
from matplotlib.path import Path as MplPath
from matplotlib.transforms import Bbox
import numpy as np
from osgeo import gdal, ogr, osr
REFERENCE = Path(__file__).resolve().parents[1] / 'figure_data'
CACHE = REFERENCE / 'bivariate_hq_display_grid.npz'
WIDTH, HEIGHT = 183, 190
LAND_COLOR = '#E5E4DE'
SLOPE_LIMIT = .003
GROUPS = ('All species', 'Threatened species')
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],'font.size':7,
    'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black',
    'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'xtick.major.width':.6,'ytick.major.width':.6,
    'xtick.major.size':2.5,'ytick.major.size':2.5,'pdf.fonttype':42,'svg.fonttype':'none','legend.frameon':False})
COLORS={'All species':'#0072B2','Threatened species':'#D55E00','Global land':'#444444'}
_M1 = np.array([[0.4122214708, 0.5363325363, 0.0514459929],
                [0.2119034982, 0.6806995451, 0.1073969566],
                [0.0883024619, 0.2817188376, 0.6299787005]])
_M2 = np.array([[0.2104542553, 0.7936177850, -0.0040720468],
                [1.9779984951, -2.4285922050, 0.4505937099],
                [0.0259040371, 0.7827717662, -0.8086757660]])
_M1I = np.linalg.inv(_M1)
_M2I = np.linalg.inv(_M2)




def _hex_to_rgb(hex_str):
    h = hex_str.lstrip('#')
    return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)], dtype=float) / 255.0

def _srgb_to_linear(c):
    c = np.asarray(c, dtype=float)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)

def _linear_to_srgb(c):
    c = np.clip(c, 0.0, 1.0)
    return np.where(c <= 0.0031308, 12.92 * c, 1.055 * np.power(c, 1 / 2.4) - 0.055)

def _rgb_to_oklab(rgb):
    lin = _srgb_to_linear(rgb)
    lms = lin @ _M1.T
    return np.cbrt(lms) @ _M2.T

def _oklab_to_linear(lab):
    lms = (lab @ _M2I.T) ** 3
    return lms @ _M1I.T

def _fit_gamut(lab):
    """Keep lightness and hue; shrink chroma until each colour is inside sRGB."""
    L = lab[..., :1]
    ab = lab[..., 1:]
    lo = np.zeros(lab.shape[:-1] + (1,))
    hi = np.ones_like(lo)
    inside = np.all((_oklab_to_linear(lab) >= -1e-6) & (_oklab_to_linear(lab) <= 1 + 1e-6),
                    axis=-1, keepdims=True)
    lo = np.where(inside, 1.0, lo)
    for _ in range(30):
        mid = (lo + hi) / 2
        test = np.concatenate([L, ab * mid], axis=-1)
        lin = _oklab_to_linear(test)
        ok = np.all((lin >= -1e-6) & (lin <= 1 + 1e-6), axis=-1, keepdims=True)
        lo = np.where(ok, mid, lo)
        hi = np.where(ok, hi, mid)
    return np.concatenate([L, ab * lo], axis=-1)

def title(ax,letter,text):
    ax.set_title(f'{letter}  {text}',loc='left',fontweight='bold',fontsize=8,pad=8)

def density(ax,x,y,identity=False):
    image=ax.hexbin(x,y,gridsize=38,mincnt=1,bins='log',cmap='Blues',linewidths=0,rasterized=True)
    if identity:
        lo=min(float(x.min()),float(y.min()));hi=max(float(x.max()),float(y.max()));pad=max((hi-lo)*.04,.00001)
        ax.set_xlim(lo-pad,hi+pad);ax.set_ylim(lo-pad,hi+pad)
        ax.plot([lo-pad,hi+pad],[lo-pad,hi+pad],color='#555555',lw=.7,ls='--',zorder=0)
    ax.xaxis.set_major_locator(MaxNLocator(4));ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'));ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    return image

def count_colorbar(fig,image,cax):
    bar=fig.colorbar(image,cax=cax,orientation='horizontal')
    maximum=float(image.get_clim()[1])
    if not np.isfinite(maximum) or maximum<=0:
        raise ValueError('A species-count color scale requires positive counts')
    highest=max(0,int(np.floor(np.log10(maximum))))
    exponents=sorted({0,highest//2,highest})
    bar.set_ticks([10**power for power in exponents])
    bar.ax.xaxis.set_major_formatter(FuncFormatter(lambda value,position:f'{value:,.0f}'))
    bar.ax.minorticks_off()
    bar.set_label('Species per hexagon',fontsize=7,labelpad=3)
    bar.outline.set_visible(False)
    return bar

def axis_mm(fig, box):
    x, y, w, h = box
    width, height = fig.get_size_inches() * 25.4
    return fig.add_axes((x / width, y / height, w / width, h / height))

def paired_panels(fig, slopes, effects, *, row=(20, 32), bar_box=(15, 8.5, 36, 1.55),
                  panel_title=title):
    pair_ax = axis_mm(fig, (15, row[0], 42, row[1]))
    data = slopes[['mean_hq_slope', 'q05_hq_slope']].to_numpy() * 1000
    if not np.isfinite(data).all():
        raise ValueError('Expected finite paired slopes for every global species')
    image = density(pair_ax, data[:, 0], data[:, 1], True)
    image.set_cmap(LinearSegmentedColormap.from_list('species_counts', ['#c9cbce', '#313237']))
    pair_ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    pair_ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    panel_title(pair_ax, 'e', 'Paired species slopes')
    pair_ax.set(xlabel='Mean HQ slope (×10⁻³ yr⁻¹)', ylabel='5th-percentile HQ slope\n(×10⁻³ yr⁻¹)')
    bar = count_colorbar(fig, image, axis_mm(fig, bar_box))
    bar.ax.tick_params(length=2, pad=1)

    realm_ax = axis_mm(fig, (86, row[0], 91, row[1]))
    realms = effects[(effects.effect == 'q05_minus_mean_paired_slope')
                     & effects.scope.str.startswith('Realm:') & (effects.scope != 'Realm:Unclassified')]
    scopes = sorted(realms.scope.unique(), key=lambda x: (x.endswith('Antarctica'), x))
    for group, shift, marker in zip(GROUPS, (-.14, .14), ('o', 's')):
        values = realms[realms.group == group].set_index('scope').loc[scopes]
        for i, row in enumerate(values.itertuples()):
            y = i + shift
            realm_ax.plot([row.ci95_low * 1000, row.ci95_high * 1000], [y, y],
                          color=COLORS[group], lw=1.08, solid_capstyle='round')
            realm_ax.plot(row.estimate * 1000, y, marker, color=COLORS[group], ms=3)
    n_antarctica = realms[(realms.scope == 'Realm:Antarctica') & (realms.group == GROUPS[0])].n_species.item()
    labels = [s.split(':', 1)[1] for s in scopes]
    labels[-1] += f' (n = {n_antarctica})'
    realm_ax.set(yticks=np.arange(len(scopes)), yticklabels=labels, ylim=(len(scopes)-.55, -.45))
    realm_ax.axvline(0, color='#969696', ls=(0, (3, 2)), lw=.62, zorder=0)
    realm_ax.tick_params(axis='y', length=0)
    realm_ax.xaxis.set_major_locator(MultipleLocator(.5))
    realm_ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    panel_title(realm_ax, 'f', 'Differences across biogeographic realms')
    realm_ax.set_xlabel('5th-percentile − mean HQ slope (×10⁻³ yr⁻¹)')
    return pair_ax, realm_ax

def points(geometry):
    return np.asarray(geometry.GetPoints())[:, :2]

def spatial_reference(code):
    ref = osr.SpatialReference()
    ref.SetFromUserInput(code)
    ref.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return ref

def open_basemap(name):
    """Use native GDAL Parquet support when available, or read identical WKB with PyArrow."""
    source = REFERENCE / (name + '.parquet')
    if ogr.GetDriverByName('Parquet') is not None:
        return ogr.Open(str(source))
    import pyarrow.parquet as pq
    table = pq.read_table(source)
    dataset = ogr.GetDriverByName('Memory').CreateDataSource('')
    layer = dataset.CreateLayer(name, srs=spatial_reference('EPSG:8857'), geom_type=ogr.wkbUnknown)
    if 'bdytyp' in table.column_names:
        layer.CreateField(ogr.FieldDefn('bdytyp', ogr.OFTInteger))
    for row in table.to_pylist():
        feature = ogr.Feature(layer.GetLayerDefn())
        geometry = ogr.CreateGeometryFromWkb(row['geometry'])
        if name == 'un_boundaries_equal_earth':
            geometry = ogr.ForceToMultiLineString(geometry)
        feature.SetGeometry(geometry)
        if row.get('bdytyp') is not None:
            feature.SetField('bdytyp', int(row['bdytyp']))
        layer.CreateFeature(feature)
    return dataset


def footprint():
    data = open_basemap('earth_footprint')
    feature = data.GetLayer(0).GetNextFeature()
    outline = points(feature.GetGeometryRef().GetGeometryRef(0))
    return outline, (*outline.min(axis=0), *outline.max(axis=0))

def geography(ax, outline, key_bounds):
    """Use the reference figure's separate context layers and UN line styles."""
    transform = osr.CoordinateTransformation(spatial_reference('EPSG:4326'),
                                             spatial_reference('EPSG:8857'))
    grid_lines = []
    for lon in range(-150, 180, 30):
        y = np.linspace(-90, 90, 361)
        grid_lines.append(np.asarray(transform.TransformPoints(list(zip(np.full_like(y, lon), y))))[:, :2])
    for lat in range(-60, 90, 30):
        x = np.linspace(-180, 180, 721)
        grid_lines.append(np.asarray(transform.TransformPoints(list(zip(x, np.full_like(x, lat)))))[:, :2])
    if key_bounds is not None:
        xmin, ymin, xmax, ymax = key_bounds
        for line in grid_lines:
            x, y = line.T
            line[(x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)] = np.nan
    ax.add_collection(LineCollection(grid_lines, colors='#ccd1d5', linewidths=.22, alpha=.65, zorder=2))

    borders = open_basemap('un_boundaries_equal_earth')
    groups = {}
    for feature in borders.GetLayer(0):
        style = feature.GetField('bdytyp')
        geometry = feature.GetGeometryRef()
        groups.setdefault(style, []).extend(points(part) for part in geometry if part.GetPointCount() >= 2)
    for style, segments in groups.items():
        line = ':' if style in (2, 4, 99) else '--' if style == 3 else '-'
        ax.add_collection(LineCollection(segments, colors='#858489' if style else '#a29da0',
                                         linewidths=.18 if style else .20,
                                         linestyles=line, alpha=.69, zorder=3))
    ax.add_patch(Polygon(outline, fill=False, ec='#929ca7', lw=.37, zorder=4))

def small_title(ax, letter, heading):
    labels = {'b': 'Mean HQ', 'c': '5th-percentile HQ', 'd': 'Low-quality habitat',
              'e': 'Paired species slopes', 'f': 'Biogeographic realms'}
    ax.set_title(labels[letter], loc='left', x=.095 if letter != 'f' else .045,
                 weight='semibold', fontsize=7.5, pad=6)
    ax.annotate(letter, (0, 1), (0, 4.7), xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline', fontsize=9, weight='bold', annotation_clip=False)

def projected_colors(grid, bounds, lut, land_color=None):
    counts, quality, slope = (grid[name] for name in ('counts', 'baseline', 'slope'))
    valid = counts > 0
    if not np.all(np.isfinite(slope[valid]) & (quality[valid] >= 0) & (quality[valid] <= 1)):
        raise ValueError('Unexpected values on modeled land')
    weights = counts[valid]
    abs_slopes = np.abs(slope[valid])
    order = np.argsort(abs_slopes)
    area = int(np.sum(weights, dtype=np.int64))
    p95 = abs_slopes[order[np.searchsorted(np.cumsum(weights[order], dtype=np.int64), .95 * area)]]
    if np.ceil(p95 * 1000) / 1000 != SLOPE_LIMIT:
        raise ValueError('Recheck displayed endpoint against area-weighted 95th percentile')
    stats = dict(land_pixels_1km=area, weighted_p95_absolute_slope=float(p95),
                 limit=SLOPE_LIMIT, lower_endpoint_percent=100 * int(counts[slope < -SLOPE_LIMIT].sum()) / area,
                 upper_endpoint_percent=100 * int(counts[slope > SLOPE_LIMIT].sum()) / area)
    west, east, south, north = grid['extent']
    nr, nc = counts.shape
    driver = gdal.GetDriverByName('MEM')
    source = driver.Create('', nc, nr, 3, gdal.GDT_UInt16)
    source.SetGeoTransform((west, (east - west) / nc, 0, north, 0, (south - north) / nr))
    source.SetProjection(spatial_reference('ESRI:54009').ExportToWkt())
    source.GetRasterBand(1).WriteArray(np.rint(quality * 256).astype('uint16'))
    trend = np.clip(slope / SLOPE_LIMIT, -1, 1)
    source.GetRasterBand(2).WriteArray(np.rint((trend + 1) * 128).astype('uint16'))
    source.GetRasterBand(3).WriteArray(counts)
    xmin, ymin, xmax, ymax = bounds
    width = 4800
    height = round(width * (ymax - ymin) / (xmax - xmin))
    display = gdal.Warp('', source, format='MEM', dstSRS='EPSG:8857', outputBounds=bounds,
                        width=width, height=height, resampleAlg='near', errorThreshold=0,
                        multithread=True, warpOptions=['NUM_THREADS=4', 'INIT_DEST=0'])
    u, v, displayed_counts = display.ReadAsArray()
    basemap = open_basemap('un_land_equal_earth')
    land = driver.Create('', width, height, 1, gdal.GDT_Byte)
    land.SetGeoTransform(display.GetGeoTransform())
    land.SetProjection(display.GetProjection())
    gdal.RasterizeLayer(land, [1], basemap.GetLayer(0), burn_values=[1])
    context = land.ReadAsArray().astype(bool)
    support = displayed_counts > 0
    rgb = np.full((height, width, 3), 255, dtype=np.uint8)
    rgb[context] = np.rint(255 * np.array(to_rgb(land_color or LAND_COLOR))).astype('uint8')
    rgb[support] = lut[v[support], u[support]]
    return rgb, support | context, stats

def position_key(fig, world, mask, outline, box, lut):
    """Check the actual swatch and every text extent against islands and coastline."""
    x, y, w, h = box
    xmin, xmax = world.get_xlim()
    ymin, ymax = world.get_ylim()
    projected = (outline - (xmin, ymin)) / (xmax - xmin, ymax - ymin) * (w, h) + (x, y)
    frame = MplPath(projected, closed=True)
    size = 15.0
    k = axis_mm(fig, (0, 0, size, size))
    k.imshow(lut, extent=(0, 1, -3, 3), origin='lower', aspect='auto', interpolation='nearest')
    k.set(xticks=[0, .5, 1], xticklabels=['0', '0.5', '1'], yticks=[-3, 0, 3],
          yticklabels=['−3', '0', '+3'], xlabel='HQ in 2000', ylabel='Trend (0.001/yr)')
    k.tick_params(length=1.8, pad=1.0)
    k.xaxis.labelpad = 1.2
    k.yaxis.labelpad = 1.0
    for spine in k.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(.35)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    mm = 25.4 / fig.dpi
    texts = [k.xaxis.label, k.yaxis.label, *k.get_xticklabels(), *k.get_yticklabels()]
    parts = [k.get_window_extent().padded(1), *(t.get_window_extent(renderer) for t in texts)]
    local = [np.array(b.extents) * mm for b in parts]
    nr, nc = mask.shape
    # Integral image for quick, exact tests at raster resolution.
    summed = np.zeros((nr + 1, nc + 1), dtype=np.int64)
    np.cumsum(mask, axis=0, dtype=np.int64, out=summed[1:, 1:])
    np.cumsum(summed[1:, 1:], axis=1, out=summed[1:, 1:])
    def clear(rect):
        xa, ya, xb, yb = rect.extents
        c0, c1 = np.clip([np.floor((xa-x)*nc/w), np.ceil((xb-x)*nc/w)], 0, nc).astype(int)
        r0, r1 = np.clip([np.floor((y+h-yb)*nr/h), np.ceil((y+h-ya)*nr/h)], 0, nr).astype(int)
        occupied = summed[r1, c1] - summed[r0, c1] - summed[r1, c0] + summed[r0, c0]
        return not occupied and not frame.intersects_bbox(rect, filled=False)
    # Prefer a key inside the southern Pacific, close to the lower-left tangent.
    candidates = [(i, j) for i in np.arange(x+15, x+42, .5) for j in np.arange(y+3, y+15, .5)]
    candidates.sort(key=lambda p: abs(p[0]-(x+31)) + 1.4*abs(p[1]-(y+8)))
    for left, bottom in candidates:
        rects = [Bbox.from_extents(*(rect + (left, bottom, left, bottom))).padded(.4) for rect in local]
        if frame.contains_point((left + size/2, bottom + size/2)) and all(clear(r) for r in rects):
            k.set_position([left / WIDTH, bottom / HEIGHT, size / WIDTH, size / HEIGHT])
            full = Bbox.union(rects)
            return full.extents, (left, bottom, size, size)
    raise ValueError('No unobstructed lower-left map legend location')

def new_figure():
    family = font_manager.FontProperties(fname=font_manager.findfont('Arial')).get_name()
    plt.rcParams.update({'font.family': family, 'font.sans-serif': [family], 'font.size': 7,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.linewidth': .46,
                         'xtick.major.width': .45, 'ytick.major.width': .45, 'xtick.major.size': 2,
                         'ytick.major.size': 2, 'axes.labelpad': 2.7})
    return plt.figure(figsize=(WIDTH / 25.4, HEIGHT / 25.4), facecolor='white')

def map_panel(fig, lut, land_color, *, baseline_table=False):
    outline, bounds = footprint()
    xmin, ymin, xmax, ymax = bounds
    with np.load(CACHE) as grid:
        rgb, occupied, stats = projected_colors(grid, bounds, lut, land_color)
    annual = pd.read_csv(REFERENCE.parent / 'tables/background_annual.csv')
    if stats['land_pixels_1km'] != annual.query("scope == 'Global' and year == 2000").grid_cells.item():
        raise ValueError('Unexpected modeled land area')
    w = 174
    box = ((WIDTH - w) / 2, 99.5, w, w * (ymax - ymin)/(xmax - xmin))
    a = axis_mm(fig, box)
    a.imshow(rgb, extent=(xmin, xmax, ymin, ymax), origin='upper', interpolation='nearest', zorder=1)
    a.images[0].set_clip_path(Polygon(outline, transform=a.transData))
    a.set(xlim=(xmin,xmax), ylim=(ymin,ymax), aspect='equal')
    a.set_axis_off()
    key_rect, key_box = position_key(fig, a, occupied, outline, box, lut)
    xa, ya, xb, yb = key_rect
    if baseline_table:
        # Reserve clear ocean directly above the key, including the values column.
        kx, ky, _, side = key_box
        nr, nc = occupied.shape
        candidates = sorted(((dx, dy) for dx in np.arange(-12, 8, .5)
                             for dy in np.arange(-.5, 8, .5)),
                            key=lambda shift: abs(shift[0]) + 1.5 * abs(shift[1]))
        for dx, dy in candidates:
            left, top = kx - 5.8 + dx, ky + side + 13.8 + dy
            l, r, b, t = left - .8, left + 22, top - 12, top + 1.8
            rows = slice(int((box[1]+box[3]-t)*nr/box[3]), int(np.ceil((box[1]+box[3]-b)*nr/box[3])))
            cols = slice(int((l-box[0])*nc/w), int(np.ceil((r-box[0])*nc/w)))
            if not occupied[rows, cols].any():
                stats['baseline_anchor_mm'] = [left, top]
                xa, xb, yb = min(xa,l), max(xb,r), max(yb,t)
                break
        else:
            raise ValueError('No clear ocean above the bivariate key for the baseline values')
    avoid = (xmin+(xa-box[0])/w*(xmax-xmin), ymin+(ya-box[1])/box[3]*(ymax-ymin),
             xmin+(xb-box[0])/w*(xmax-xmin), ymin+(yb-box[1])/box[3]*(ymax-ymin))
    geography(a, outline, avoid)
    for x, label, size in [(4.5,'a',9),(9,'Habitat quality in 2000 and its 2000–2020 trend',7.5)]:
        fig.text(x/WIDTH,187/HEIGHT,label,fontsize=size,weight='bold',va='center')
    stats.update(map_box_mm=box, legend_swatch_mm=key_box, legend_full_mm=key_rect.tolist())
    return a, stats



