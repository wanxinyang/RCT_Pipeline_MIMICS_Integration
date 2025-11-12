# RCT2MIMICS: Integrating LiDAR-derived Tree Structure with the Michigan Microwave Scattering Model

A modular, three-step workflow that integrates 3D tree structural data from the [RCT-pipeline](https://github.com/wanxinyang/rct-pipeline) with the [Michigan Microwave Scattering Model (MIMICS)](https://codeocean.com/capsule/8976769/tree/v1) for simulating radar backscatter from LiDAR-scanned trees.


## Overview

**RCT2MIMICS** provides an automated workflow to bridge LiDAR-derived structural reconstructions and radar scattering simulations.  

The process comprises three main stages:

1. **Extract tree-level attributes** from RCT-derived outputs.  
2. **Generate MIMICS input parameters** from reconstructed tree and branch structural attributes.  
3. **Run the MIMICS model** using generated parameters, enabling parallelised simulations and automated results parsing.

```
Step 1: extract_rct_tree_attrs.py
    ↓ (tree attributes CSV)
Step 2: create_parameters.py  
    ↓ (MIMICS input parameters CSV)
Step 3: run_model.py
    ↓ (MIMICS backscatter results CSV)
```

---

## Prerequisites

1. **TLS Data Registration**: 
    - Scan position matrices (`matrix/*.DAT`)

2. **RCT Pipeline**: Process your TLS data using [RCT-pipeline](https://github.com/wanxinyang/rct-pipeline) to generate individual tree segment files with the expected directory structure as follows:
    ```
    rct_extraction/
    ├── {country}_{plot}_raycloud_{tile_x}_{tile_y}_trees_info.txt
    ├── {country}_{plot}_raycloud_{tile_x}_{tile_y}_trees.txt  
    ├── {country}_{plot}_raycloud_{tile_x}_{tile_y}_treesplit/
    │   ├── {country}_{plot}_raycloud_{tile_x}_{tile_y}_trees_{tree_id}_info.txt
    │   ├── {country}_{plot}_raycloud_{tile_x}_{tile_y}_segmented_{tree_id}_leaves.ply
    │   └── ...
    └── tile_index.dat
    ```


3. **MIMICS Model**: Compile MIMICS 
    ```bash
    # make the script executable
    $ chmod +x model/code/run.sh
    # compile and build mimics1.5
    $ ./run.sh
    ```

---

## Step 1: Extract Tree Attributes

**Script**: `extract_rct_tree_attrs.py`

### Example Usage

```bash
python extract_rct_tree_attrs.py \
    -c [COUNTRY_NAME] \
    -p [PLOT_ID] \
    -d [YYYY-MM-DD] \
    -r /path/to/rct_extraction/ \
    -m /path/to/matrix/ \
    -t /path/to/tile_index.dat \
    -o /path/to/results/ \
    -f /path/to/figures/ 
```

### Parameters Explanation

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-c, --country` | Country code (e.g., angola, gabon) | Required |
| `-p, --plotid` | Plot ID (e.g., P02, P09) | Required |
| `-d, --date` | Survey date (YYYY-MM-DD) | None |
| `-r, --dir_rct_extraction` | Directory containing RCT extraction results | Required |
| `-m, --dir_matrix` | Directory with scan position matrices (*.DAT) | Required |
| `-t, --fp_tile` | RCT tile index file path | Required |
| `-o, --odir` | Output directory for CSV files | Current directory |
| `-f, --figdir` | Output directory for figures | Current directory |
| `--rct_tile_len` | RCT tile length in meters | 50 |
| `--rct_tile_overlap` | RCT tile overlap in meters | 10 |
| `--n_cores` | Number of parallel workers for leaf processing | CPU count - 1 |

### Output

```
{date}_{country}_{plotid}_treeLevel_attr_rct_tile{rct_tile_len}m_overlap{rct_tile_overlap}m.csv
```

**Columns include**: `plot_id`, `tile`, `tree_id`, `x`, `y`, `z`, `height`, `DBH`, `crown_radius`, `in_plot`, `leaf_area_m2`, `num_leaves`, and more.



---

## Step 2: Generate MIMICS Parameters

**Script**: `create_parameters.py`


**Note**: This script can work with or without census data:
- **With census matching** (`-c` provided): Processes only census-matched trees with species information
- **Without census matching** (`-c` not provided): Processes all in-plot trees from tree attributes CSV

### Example Usage

```bash
python create_parameters.py \
    -t /path/to/tree_attributes.csv \
    -r /path/to/rct_extraction/ \
    -c /path/to/census_matched_stems.csv \
    -w /path/to/wood_density.csv \
    -odir /path/to/results/ \
    --min_branch_order 3 \
    --max_branch_order 4 \
    -d 10.0 \
    -H 3.0 \
    -n 8 \
    -f 0.43 1.2 \
    -s 0.25 0.5 0.75 \
    --canopy_density 0.015 0.06 0.24 0.48 0.72
```

### Parameters Explanation

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-t, --tree_attr_csv` | Tree attributes CSV from Step 1 | Required |
| `-r, --rct_extraction_dir` | Directory with RCT segment files | Required |
| `-c, --census_matched_stems_csv` | Census to RCT mapping CSV | Optional |
| `-w, --wood_density_csv` | Wood density lookup CSV | Optional |
| `-odir, --output_dir` | Output directory for parameter CSV | Current directory |
| `-a, --area` | Plot area in hectares (for density calculation) | 1.0 |
| `--min_branch_order` | Minimum branch order required | 3 |
| `--max_branch_order` | Maximum branch order to process | 4 |
| `-d, --DBH_threshold` | Minimum DBH threshold (cm) | 10.0 |
| `-H, --H_threshold` | Minimum height threshold (m) | 3.0 |
| `-n, --n_cores` | Number of parallel workers | CPU count - 1 |
| `-f, --radar_frequency` | Radar frequencies for sweep (GHz, space-separated) | 0.43 1.2 5.4 |
| `--canopy_density` | Canopy density values (trees/m², space-separated) | Auto-calculated |
| `-s, --soil_moisture` | Soil moisture values for sweep (space-separated) | 0.25 0.5 0.75 |

#### Parameter Sensitivity Sweeps

The script generates parameter combinations across multiple values for sensitivity analysis:

```python
# Default sweep ranges (can be customised via command-line arguments)
-f, --radar_frequency:     [0.43, 1.2, 5.4]         
--canopy_density:          Auto-calculated or custom 
-s, --soil_moisture:       [0.25, 0.5, 0.75]        

# Example: 3 frequencies × 5 densities × 3 soil moisture = 45 parameter sets per tree
```

#### Data-Derived Parameters

Automatically calculated from RCT data:
- **Branch statistics** (per order 1-4): length, diameter, volume, surface area, density
- **Branch orientation PDFs**: Best-fit MIMICS PDF type using KL divergence (18 types available)
- **Trunk properties**: length, diameter, surface area, volume
- **Crown properties**: thickness, radius, volume
- **Leaf properties**: total area, leaf density

#### Fixed Parameters

Set in the script (can be customised):

```python
# Sensor parameters
angle = 30  # degrees - radar incidence angle

# Ground parameters  
rms_height = 1.5  # cm - surface roughness
correlation_length = 17.5  # cm
percent_sand = 40  # %
percent_clay = 10  # %

# Moisture content (same across all branch orders)
trunk_moisture = 0.5  # gravimetric fraction
branch_moisture = 0.5  # applied to all branch orders

# Wood density (species-specific from wood_density.csv or default)
default_wood_density = 0.5  # g/cm³
```

### Output

```
mimics_inputs_{date}_{country}_{plotid}.csv
```

Each row represents one complete MIMICS parameter set. Columns include:
- Tree metadata (species, plot_id, tile_id, tree_id)
- Sensor parameters (frequency, angle)
- Ground parameters (soil_moisture, roughness, composition)
- Tree structure (trunk dimensions, crown properties)
- Branch parameters (length, diameter, density, PDF type for orders 1-5)
- Leaf parameters (density)

---

## Step 3: Execute MIMICS Model

**Script**: `run_model.py`

### Example Usage

```bash
python run_model.py \
    -i /path/to/mimics_inputs.csv \
    -c model/data/ \
    -s model/code/ \
```

### Parameters Explanation

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-i, --input` | MIMICS input parameters CSV from Step 2 | Required |
| `-c, --mimics-config` | MIMICS configuration/data directory | `./model/data` |
| `-s, --mimics-src` | MIMICS source directory (with executable) | `./model/code` |
| `-o, --output` | Output CSV file path | Auto-generated* |
| `--n_cores` | Number of parallel workers | CPU count - 2 |
| `--retries` | Retry attempts for failed runs | 3 |
| `--preserve-individual-results` | Save individual tree results | False |
| `--debug` | Enable debug-level logging | False |

\* If not specified, replaces `mimics_inputs_` with `mimics_outputs_` in input filename



### Output
**Main results**
```
mimics_outputs_{date}_{country}_{plotid}.csv
```

Contains all input parameters plus MIMICS results:
- Backscatter components (σ₀): total, HH, VV, HV
- Mechanism breakdown: ground, trunk, canopy contributions

**Individual results** (if `--preserve_individual_results` enabled):
```
individual_results/{country}_{date}_{plotid}_{tile}_Tree_{tree_id}/
```

---

## Citation

If you use **RCT2MIMICS** or its script(s) in your research or publications, please cite this repository:

> Yang, W. and Thomas, R. (2025) ‘RCT2MIMICS: Integrating LiDAR-derived tree structure with the Michigan Microwave Scattering Model’. Zenodo. doi:10.5281/zenodo.17594216.

**BibTeX**
```bibtex
@misc{yang2025rct2mimics,
  author       = {Wanxin Yang and Russell Thomas},
  title        = {RCT2MIMICS: Integrating LiDAR-derived tree structure with the Michigan Microwave Scattering Model},
  year         = {2025},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.17594216},
  url          = {https://zenodo.org/record/17594216}
}
```

---

## To-Do

- [ ] Expose Step 2 fixed parameters (angle, moisture, ground properties) as CLI arguments
- [ ] Create Step 1.5: Tree selection script to filter trees before parameter generation (separate filtering logic from `create_parameters.py`)

---
