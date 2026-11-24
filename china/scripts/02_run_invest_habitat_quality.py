"""Run the InVEST Habitat Quality model for multiple years.

The script builds one annual threats table per year using threat raster paths
provided by filename templates. It avoids local absolute paths and writes all
model outputs under a user-specified output directory.
"""

from __future__ import annotations

import argparse
import csv
import importlib
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


def write_threats_table(year: int, threat_dir: Path, output_csv: Path) -> Path:
    """Create an annual InVEST threats table."""
    threat_paths = {
        "urb": threat_dir / f"urban_{year}.tif",
        "crp": threat_dir / f"cropland_{year}.tif",
        "rr": threat_dir / f"roads_{year}.tif",
    }
    missing = [str(path) for path in threat_paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing threat rasters:\n" + "\n".join(missing))

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as file_obj:
        writer = csv.writer(file_obj)
        writer.writerow(["THREAT", "MAX_DIST", "WEIGHT", "DECAY", "CUR_PATH"])
        with (Path(__file__).resolve().parents[1] / "parameters/threats.csv").open(encoding="utf-8-sig") as source:
            parameters = list(csv.DictReader(source))
        for params in parameters:
            threat_name = params["THREAT"]
            writer.writerow([
                threat_name,
                float(params["MAX_DIST"]) / 1000,  # archived metres to InVEST 3.14.2 kilometres
                params["WEIGHT"],
                params["DECAY"],
                threat_paths[threat_name].resolve(),
            ])
    return output_csv


def build_invest_args(
    year: int,
    lulc_dir: Path,
    threat_dir: Path,
    sensitivity_table: Path,
    output_dir: Path,
    half_saturation: float,
) -> dict[str, object]:
    """Build the InVEST Habitat Quality argument dictionary for one year."""
    workspace = output_dir / f"HQ_{year}"
    workspace.mkdir(parents=True, exist_ok=True)
    lulc_path = lulc_dir / f"land_cover_{year}.tif"
    if not lulc_path.exists():
        raise FileNotFoundError(f"Missing land-cover raster: {lulc_path}")
    threats_csv = write_threats_table(year, threat_dir, workspace / f"threats_{year}.csv")
    # InVEST 3.14.2 uses lulc as its sensitivity-table key.
    local_sensitivity = workspace / "sensitivity_invest3142.csv"
    with sensitivity_table.open(encoding="utf-8-sig") as source, local_sensitivity.open("w", newline="", encoding="utf8") as dest:
        reader = csv.reader(source)
        writer = csv.writer(dest)
        writer.writerow(["lulc" if c.lower() == "lucode" else c.lower() for c in next(reader)])
        writer.writerows(reader)
    return {
        "workspace_dir": str(workspace),
        "results_suffix": str(year),
        "lulc_cur_path": str(lulc_path),
        "threats_table_path": str(threats_csv),
        "sensitivity_table_path": str(local_sensitivity.resolve()),
        "half_saturation_constant": half_saturation,
        "n_workers": -1,
    }


def run_year(params: tuple[int, Path, Path, Path, Path, float]) -> str:
    """Run one annual InVEST model task."""
    year, lulc_dir, threat_dir, sensitivity_table, output_dir, half_saturation = params
    args = build_invest_args(year, lulc_dir, threat_dir, sensitivity_table, output_dir, half_saturation)
    habitat_quality = importlib.import_module("natcap.invest.habitat_quality")
    habitat_quality.execute(args)
    return f"Completed {year}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lulc-dir", type=Path, required=True, help="Directory with land_cover_YEAR.tif rasters.")
    parser.add_argument("--threat-dir", type=Path, required=True, help="Directory with urban_YEAR.tif, cropland_YEAR.tif, roads_YEAR.tif.")
    parser.add_argument("--sensitivity-table", type=Path, default=Path(__file__).resolve().parents[1] / "parameters/sensitivity.csv")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--start-year", type=int, default=2000)
    parser.add_argument("--end-year", type=int, default=2020)
    parser.add_argument("--half-saturation", type=float, default=0.5)
    parser.add_argument("--workers", type=int, default=2)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    tasks = [
        (year, args.lulc_dir, args.threat_dir, args.sensitivity_table, args.output_dir, args.half_saturation)
        for year in range(args.start_year, args.end_year + 1)
    ]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(run_year, task) for task in tasks]
        for future in as_completed(futures):
            print(future.result())


if __name__ == "__main__":
    main()
