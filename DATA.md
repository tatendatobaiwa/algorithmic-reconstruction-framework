# Data sources and licences

| Data | Source | Used for | In this repository |
|---|---|---|---|
| Orbital element sets (TLE/OMM) | Space-Track, 18th and 19th Space Defense Squadrons | Orbit propagation of every satellite for each window | No. Not redistributed under the Space-Track user agreement; each run records the epochs it used in its manifest. |
| Protected-area boundaries | UNEP-WCMC and IUCN, Protected Planet: The World Database on Protected Areas (WDPA) | Region polygons, discretised to H3 resolution 7 | No. The WDPA terms allow non-commercial use and do not allow redistribution. Files listing the grid cells inside the regions are held back until UNEP-WCMC confirms that they may be published. |
| Sentinel-2 L2A scene classification | Copernicus programme, read from the public cloud archive | Observed cloud in amendment A8 | Derived per-overpass statistics only (`results/amendments/a8results/`). Contains modified Copernicus Sentinel data 2026. |
| Sensor parameters (swath, sun limits, data rates) | Operator documentation, cited in Cell 0 of the pipeline | Footprints and data volume | Values and sources are in the code. |

Everything under `results/` is derived by the code in this repository and is released with it under the MIT licence.
