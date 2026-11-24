# Global model inputs and parameters

The global extension includes China. The source roster uses terrestrial, extant, native or reintroduced records for resident, breeding or non-breeding ranges. The five accepted Red List categories are LC, NT, VU, EN and CR. The analytical sample requires overlap with common valid land. The species ID links the source roster, status table and annual observations.

The production model uses InVEST 3.14.2, a 1,000 m Mollweide grid (ESRI:54009), six-degree processing tiles and a 15 km halo. Built-up land, cropland and roads have maximum influence distances of 10, 2 and 3 km, respectively, weights of 1, 0.25 and 1, and exponential, linear and linear decay. The half-saturation constant is 0.5. The code converts the configured distances in metres to the units required by the model. The China threat table retains its original metre values.

`original_parameters.csv` gives the originating habitat scores and sensitivities. The ESA crosswalk lists each annual land-cover code, its parameter source or mixture, and the two land-cover threat flags. `sensitivity_invest3142.csv` supplies the model-ready categories. ESA CCI v2.0.7cds is used through 2015 and C3S v2.1.1 thereafter. The NetCDF metadata records the inspected 2000, 2010 and 2020 files. GRIP4 roads are static. The global model has no China-style railway layer.

## Running the upstream model

Use a separate environment with InVEST 3.14.2, GDAL with NetCDF, FileGDB reading and SQLite support, shapely 2, numpy, pandas, scipy and portalocker. Supply the input products from manuscript Table 1 under a directory of your choice. Set `HQ_INPUTS` to that directory. Relative source paths in `run_config.json` resolve against it. An absolute source path in your own config is also accepted. Set `HQ_CONFIG` to a local config copy and `HQ_OUTPUT` to a separate, writable run directory. Run `python -B global/scripts/run_global.py` from the package root.

| Config field | Required input |
|---|---|
| esa_archives | Directory of annual ESA ZIP archives with a single matching NetCDF member per year. An optional netcdf map may specify pre-extracted files by year. |
| roads | GRIP4 road line dataset in geographic coordinates. |
| species_root | Authorized IUCN release, grouped in subdirectories with shapefiles, corresponding to the roster in global_species_screened.csv. |
| land_boundaries | Geographic country polygons with the NAME_0 field. Exterior rings define the land footprint. |
| ecoregions | Geographic RESOLVE Ecoregions 2017 polygons with the REALM field. |
| hfp_root | Directory containing hfp2000.tif and hfp2020.tif for the auxiliary paired-endpoint comparison. |

The release uses disk scratch by default. Workers, scratch location (`HQ_SCRATCH` with the memory backend), disk limits and the random seed are configurable. Use a fresh output directory when changing inputs or parameters. The pipeline resumes completed stages within a run directory. It checks the 15 and 25 km halos on a reference tile. `model_validation_reference.json` records the numerical checks from the original production environment. Group confidence intervals use 2,000 whole-species trajectory resamples with base seed 20261117.

China upstream entry points use pre-screened, China-clipped species records. The richness script accepts the analysis raster as a template to preserve its alignment. The annual InVEST entry point expects land_cover_YEAR.tif and urban_YEAR.tif, cropland_YEAR.tif and roads_YEAR.tif. The China roads layer includes roads and railways. Supply aligned rasters, the template and authorized range files separately. The exposure entry point computes mean HQ, within-range P5 and the percentage of valid pixels at or below 0.239.

China raw GIS processing additionally uses GeoPandas, rasterio and tqdm. For a new China forest fit, install scikit-learn (1.2 or later), shap, joblib and tqdm alongside pandas and PyArrow. Run `python -B china/scripts/gwrf_shap_annual.py --output-dir ../review-output/refit` and then `python -B china/scripts/merge_annual_shap.py --input-dir ../review-output/refit/annual --output-dir ../review-output/refit`. The refit writes local forests to disk. The original tuning subsample was unseeded. The release exposes a seed for new fits and preserves the archived split.
