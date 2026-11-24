"""Extract habitat quality exposure metrics for species range polygons."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask
from tqdm import tqdm


def load_species_ranges(species_dir: Path) -> gpd.GeoDataFrame:
    """Load and concatenate all species range shapefiles in a directory."""
    shapefiles = sorted(species_dir.glob("*.shp"))
    if not shapefiles:
        raise FileNotFoundError(f"No shapefiles found in {species_dir}")
    frames = [gpd.read_file(path).to_crs('EPSG:4326') for path in shapefiles]
    ranges = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), crs='EPSG:4326')
    if ranges.id_no.isna().any():
        raise ValueError('Every range record requires a species ID')
    return ranges.dissolve(by='id_no', as_index=False)


def extract_values(geometry, raster_path: Path) -> np.ndarray:
    """Extract valid raster values under a single geometry."""
    with rasterio.open(raster_path) as src:
        try:
            out_image, _ = mask(src, [geometry], crop=True, filled=False, all_touched=False)
        except ValueError as exc:
            if "do not overlap" in str(exc):
                return np.empty(0)
            raise
        values = out_image[0].compressed()
        valid = np.isfinite(values)
        if src.nodata is not None:
            valid &= values != src.nodata
        valid &= np.abs(values) < 1e10
        return values[valid]


def exposure_metrics(values: np.ndarray, low_threshold: float) -> dict[str, float | int]:
    """Compute exposure metrics for one species-year record."""
    if values.size == 0:
        return {"q_mean": np.nan, "q_std": np.nan, "r_low": np.nan, "q_5pct": np.nan, "n_pixels": 0}
    r_low = np.nan if np.isnan(low_threshold) else float((values <= low_threshold).mean() * 100.0)
    return {
        "q_mean": float(np.mean(values)),
        "q_std": float(np.std(values, ddof=1)) if values.size > 1 else 0.0,
        "r_low": r_low,
        "q_5pct": float(np.percentile(values, 5)),
        "n_pixels": int(values.size),
    }


def process_record(args: tuple[dict, Path, float, int]) -> dict[str, object]:
    """Process one species record for one year."""
    record, raster_path, low_threshold, year = args
    values = extract_values(record["geometry"], raster_path)
    metrics = exposure_metrics(values, low_threshold)
    result: dict[str, object] = {"year": year, **metrics}
    for column in ["id_no", "binomial", "sci_name", "category", "taxon", "presence", "origin"]:
        if column in record:
            result[column] = record[column]
    return result


def process_year(year: int, species_gdf: gpd.GeoDataFrame, raster_dir: Path, workers: int) -> list[dict[str, object]]:
    """Extract exposure metrics for one year."""
    raster_path = raster_dir / f"HQ{year}.tif"
    if not raster_path.exists():
        raster_path = raster_dir / f"quality_out_{year}.tif"
    if not raster_path.exists():
        raster_path = raster_dir / f'HQ_{year}/quality_out_{year}.tif'
    if not raster_path.exists():
        raise FileNotFoundError(f'Missing habitat-quality raster for {year}: {raster_path}')

    with rasterio.open(raster_path) as src:
        target_crs = src.crs
    species_projected = species_gdf.to_crs(target_crs) if species_gdf.crs != target_crs else species_gdf
    low_threshold = pd.read_csv(Path(__file__).resolve().parents[1] / 'parameters/low_quality_threshold.csv').upper_bound.item()
    records = [(row.to_dict(), raster_path, low_threshold, year) for _, row in species_projected.iterrows()]

    results: list[dict[str, object]] = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(process_record, record) for record in records]
        for future in tqdm(as_completed(futures), total=len(futures), desc=f"Year {year}"):
            results.append(future.result())
    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-dir", type=Path, required=True)
    parser.add_argument("--raster-dir", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--start-year", type=int, default=2000)
    parser.add_argument("--end-year", type=int, default=2020)
    parser.add_argument("--workers", type=int, default=3)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    species_gdf = load_species_ranges(args.species_dir)
    all_results: list[dict[str, object]] = []
    for year in range(args.start_year, args.end_year + 1):
        all_results.extend(process_year(year, species_gdf, args.raster_dir, args.workers))
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(all_results).to_csv(args.output_csv, index=False)
    print(f"Saved exposure metrics to {args.output_csv}")


if __name__ == "__main__":
    main()
