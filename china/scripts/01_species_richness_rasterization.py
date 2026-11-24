"""Count distinct species per reference cell from authorized, China-clipped range records."""
import argparse
from pathlib import Path
import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio
from rasterio.features import rasterize


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True, help='Pre-screened China range files for the selected group.')
    parser.add_argument('--template', type=Path, required=True, help='Study raster specifying CRS, extent and cell alignment.')
    parser.add_argument('--output-raster', type=Path, required=True)
    args = parser.parse_args()
    with rasterio.open(args.template) as grid:
        frames = [gpd.read_file(f).to_crs(grid.crs)[['id_no','geometry']] for f in sorted(args.input_dir.glob('*.shp'))]
        if not frames:
            raise ValueError('No range files supplied')
        records = gpd.GeoDataFrame(pd.concat(frames,ignore_index=True),crs=grid.crs)
        if records.id_no.isna().any():
            raise ValueError('Every source range requires a species ID')
        species = records.dissolve(by='id_no')
        # Burn each dissolved species once. Separately rasterizing each species avoids
        # GDAL additive burning of a multipolygon boundary within the same cell.
        profile = dict(driver='GTiff',width=grid.width,height=grid.height,transform=grid.transform,
                       crs=grid.crs,count=1,dtype='uint16',compress='deflate',tiled=True)
        args.output_raster.parent.mkdir(parents=True,exist_ok=True)
        from rasterio.features import geometry_window
        from rasterio.errors import WindowError
        counts = np.zeros(grid.shape,dtype=np.uint16)
        for geom in species.geometry:
            if geom is None or geom.is_empty:
                continue
            try:
                window = geometry_window(grid,[geom])
            except WindowError:
                continue
            rows, cols = window.toslices()
            counts[rows,cols] += rasterize([geom],out_shape=(int(window.height),int(window.width)),
                                          transform=grid.window_transform(window),default_value=1,
                                          fill=0,dtype='uint8',all_touched=False)
        with rasterio.open(args.output_raster,'w',**profile) as output:
            output.write(counts,1)


if __name__ == '__main__':
    main()
