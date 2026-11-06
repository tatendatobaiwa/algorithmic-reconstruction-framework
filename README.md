# Algorithmic Reconstruction Framework (ARF)

Code, protocol and derived results for selecting images from Earth-observation satellites already in orbit so that a mosaic of a large protected area is completed as early as the orbits allow, with as few images as possible. The motivating use cases are fire and flood mapping in Botswana. The primary region is the Central Kalahari Game Reserve with the Okavango Delta; the Lower Zambezi (Zambia) and the Karoo Nature Reserve (South Africa) test whether the results transfer to other regions.

Author: Tatenda Allen Tobaiwa, Botswana International University of Science and Technology (BIUST). MSc research.

Status: release of record **1.5.3** for Paper 1 (real Earth-observation fleets). Releases 1.5.4 and 1.5.5 are kept for the companion study on a hypothetical dense constellation; 1.5.5 is under review and has no launch cell on purpose.

## What is in the repository

| Folder | Contents |
|---|---|
| `pipeline/v1.5.3/` | Release of record: the pipeline as one file (`algorithmic_reconstruction_framework.py`, SHA-256 `14426c89…7950`), the same code as 14 notebook cells, the Colab notebook (outputs cleared), the batch runner, the parity check and its expected values, and the validation suite (109 checks). |
| `pipeline/v1.5.4/`, `pipeline/v1.5.5/` | Later releases for the dense-constellation arm, each with release notes stating exactly what changed. |
| `protocol/` | The confirmatory protocol, frozen before any protocol window was run, with amendments A1 to A13 and the status of each when written. |
| `analyses/A4_…` to `analyses/A12_A13_…` | One folder per registered amendment: the analysis script, its Colab install cell, its tests and a README. |
| `results/summary/` | Aggregated results per window and per design (JSON and CSV), and `facts.json`, the registry that every number in the manuscript is checked against. |
| `results/amendments/` | Per-amendment results (per window and per onset). |
| `results/reference_tables/` | The manuscript tables as generated; the tests regenerate them and compare. |
| `manuscript_tools/` | Scripts that draw the figures and write the tables from `results/`. Each script asserts its inputs against the stored summaries before drawing. |
| `tests/` | Repository checks that need no network, no Orekit and no licensed data. |

The manuscript text is not included until it is submitted.

## Installation

Python 3.11 to 3.13. Orekit runs through `orekit-jpype` with a bundled JDK, so no separate Java installation is needed.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python -m pytest -q tests        # 7 checks, about 10 s
```

The tests confirm that the release file is the one the protocol was frozen against, that every Python file compiles, that the manuscript tables and figures regenerate from `results/` and equal the reference copies, and that no licensed data are present.

## Reproducing the runs

The runs were made in Google Colab with outputs written to Google Drive. To repeat them:

1. Put the files of `pipeline/v1.5.3/` in `MyDrive/ARF_code/`.
2. Add two Colab secrets, `SPACETRACK_USERNAME` and `SPACETRACK_PASSWORD`, for a [Space-Track](https://www.space-track.org) account. The code reads them only from the secrets; it never stores or prints them.
3. Obtain the WDPA boundaries of the protected areas from [Protected Planet](https://www.protectedplanet.net) under its own terms (see `DATA.md`): download each site as a shapefile and put the downloads in the folder named by `ARF_BOUNDARY_DIR`. The pipeline matches the files by WDPA site ID and checks the tiled area against the WDPA area.
4. Run the parity check (`arf_colab_parity.py`, expect `5/5 runs MATCH`), then the batch runner: 8 primary, 32 secondary and 16 exploratory runs. Details and expected outputs are in `pipeline/v1.5.3/RELEASE_NOTES.txt`.
5. For an amendment, run the install cell and script in its folder under `analyses/` on the completed run folders; each README gives the command and run time.

Main settings (from `protocol/arf_confirmatory_protocol.json`, which holds the complete list): eight dry-season windows of 48 h starting 3 August to 21 September 2026; H3 resolution 7; onsets shifted over 24 h in steps of 20 min; at most 30 selection steps; target 95% of reachable cells; cloud probability 0.35 with a correlation length of 50 km and 3 h; 20 cloud seeds; random seed 42. Sun elevation limits are set per sensor in Cell 0, each with its source.

## Data that are not in the repository

* **Orbital element sets** from Space-Track are not redistributed under its user agreement. The pipeline downloads them for each window; the epochs used are recorded in each run's manifest.
* **WDPA/WDPCA boundaries** are licensed for non-commercial use without redistribution. Files that list grid cells inside the reserves (`arf_refined_events.csv`, `arf_mosaic_tiles.csv`) are held back until UNEP-WCMC confirms that they may be published.
* **Sentinel-2 cloud masks** (amendment A8) are read from the public archive at run time. Contains modified Copernicus Sentinel data 2026.

## Licence and citation

Code under the MIT licence (`LICENSE`). Please cite using `CITATION.cff`. Statement on the use of AI tools: `AI_ASSISTANCE.md`.
