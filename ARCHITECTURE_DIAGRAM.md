# Architecture Comparison

## OLD APPROACH (Version 1.0 - Russ original)

```
┌─────────────────────────────────────────────────────────────────┐
│ Step 1: get_tree_data.py                                        │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  RCT Extraction Files (_info.txt)                               │
│            ↓                                                     │
│  Read all trees from all tiles                                  │
│            ↓                                                     │
│  segment_data/angola_p02_combined.csv  ← REDUNDANT FILE         │
│  segment_data/angola_p09_combined.csv                           │
│  segment_data/angola_p10_combined.csv                           │
│  segment_data/angola_p14_combined.csv                           │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
                             ↓
┌─────────────────────────────────────────────────────────────────┐
│ Step 2: create_parameters_update.py                             │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  Read segment_data/*.csv (all at once)                          │
│            ↓                                                     │
│  Load census data                                               │
│  Load wood density data                                         │
│            ↓                                                     │
│  Process ALL trees sequentially                                 │
│   - Calculate angles                                            │
│   - Calculate PDFs                                              │
│   - Calculate branch stats                                      │
│   - Generate parameters                                         │
│            ↓                                                     │
│  model_input.csv                                                │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
                             ↓
┌─────────────────────────────────────────────────────────────────┐
│ Step 3: run_model.py (original version)                         │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  Read model_input.csv                                           │
│            ↓                                                     │
│  For each parameter set (basic parallel):                       │
│   - Configure MIMICS inputs                                     │
│   - Run MIMICS model                                            │
│   - Parse results                                               │
│            ↓                                                     │
│  mimics_results.csv                                             │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘

═══════════════════════════════════════════════════════════════════
```


## NEW APPROACH (Version 2.0 - Wanxin refactored)

```
┌─────────────────────────────────────────────────────────────────┐
│ Step 1: extract_rct_tree_attrs.py ⭐ NEW SCRIPT                 │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  RCT Raw Outputs:                                               │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ • *_trees_info.txt (tree geometry & segments)           │   │
│  │ • *_trees.txt (tree locations)                          │   │
│  │ • *_leaves.ply (leaf meshes for area calculation)       │   │
│  │ • SOP matrix files (*.DAT - scan positions)             │   │
│  │ • RCT tile index files                                  │   │
│  └──────────────────────────────────────────────────────────┘   │
│            ↓                                                     │
│  Processing Steps:                                              │
│  1. Read & merge tree info + location data                      │
│  2. Calculate plot boundaries from scan positions               │
│  3. Filter in-plot trees (spatial boundary)                     │
│  4. Calculate leaf areas (PARALLEL processing) ⚡               │
│  5. Merge tree attributes with leaf data                        │
│            ↓                                                     │
│  Output: {date}_{country}_{plot}_treeLevel_attr_rct_*.csv        │
│          ↑ No redundant segment files!                          │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
                             ↓
┌─────────────────────────────────────────────────────────────────┐
│ Step 2: create_parameters.py ⭐ REFACTORED                      │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  Input Data Sources:                                            │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ • Tree attributes CSV (from Step 1)                     │   │
│  │ • Matched stems CSV (census ↔ RCT mapping)              │   │
│  │ • Wood density CSV (species lookup)                     │   │
│  │ • RCT segment files (_info.txt) ← Read per-tree! ⚡     │   │
│  └──────────────────────────────────────────────────────────┘   │
│            ↓                                                     │
│  For each tree (PARALLEL): ⚡                                   │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ 1. Load RCT segment data for specific tree (on-demand)   │   │
│  │ 2. Calculate branch geometry & statistics                │   │
│  │ 3. Compute branch orientation angles                     │   │
│  │ 4. Fit MIMICS PDF distributions (KL divergence)         │   │
│  │ 5. Generate parameter combinations (sensitivity)         │   │
│  │ 6. Combine tree, leaf, and wood density data            │   │
│  └──────────────────────────────────────────────────────────┘   │
│            ↓                                                     │
│  Output: mimics_inputs_{date}_{country}_{plot}.csv                     │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
                             ↓
┌─────────────────────────────────────────────────────────────────┐
│ Step 3: run_model.py ⭐ IMPROVED                                │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  Input: MIMICS parameter sets CSV (from Step 2)                 │
│            ↓                                                     │
│  For each parameter set (PARALLEL): ⚡                          │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ 1. Create isolated temporary directory                   │   │
│  │ 2. Configure MIMICS input files with tree parameters    │   │
│  │ 3. Execute MIMICS radar backscatter model               │   │
│  │ 4. Parse output results (σ₀ values)                     │   │
│  │ 5. Optional: preserve individual results                │   │
│  │ 6. Retry logic for failed runs                          │   │
│  │ 7. Progress tracking with tqdm                          │   │
│  └──────────────────────────────────────────────────────────┘   │
│            ↓                                                     │
│  Output: mimics_outputs_{date}_{country}_{plot}.csv                    │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘

  ✅ Benefits of Current 3-Step Approach:
  - Clear separation of concerns (extract → parameterise → model)
  - Comprehensive tree attribute extraction (Step 1)
  - Spatial filtering and in-plot tree identification
  - Flexible I/O with CLI
  - On-demand segment loading (Step 2, low memory)
  - Full parallelisation at all compute-intensive steps
  - Better error handling and retry logic (Step 3)
  - Progress tracking and logging throughout
  - Intermediate checkpoints for debugging
  - Reusable tree attributes (Step 1 output)
  - Flexible parameter generation (Step 2 and Step 3)
  - Scalable model execution

═══════════════════════════════════════════════════════════════════
```

## DATA FLOW DETAIL 

### Step 1: Tree Attribute Extraction

```
┌─────────────────────────────────────────────────────────────────┐
│ RCT Raw Files (use Angola Plot P02 as an example)                              │
│ ───────────────────────────────────────────────────────────────  │
│ • angola_p02_raycloud_*_trees_info.txt (tree and branch info)   │
│ • angola_p02_raycloud_*_trees.txt (tree positions)              │
│ • angola_p02_raycloud_*_leaves.ply (leaf meshes)                │
│ • scan position matrices (*.DAT files)                          │
└─────────────────────────────────────────────────────────────────┘
                       ↓
┌─────────────────────────────────────────────────────────────────┐
│ Processing (extract_rct_tree_attrs.py)                          │
│ ───────────────────────────────────────────────────────────────  │
│ 1. Read tree info files → merge with locations                  │
│ 2. Remove duplicate trees across tiles                          │
│ 3. Calculate plot boundary from scan positions                  │
│ 4. Filter trees within plot boundaries                          │
│ 5. Calculate leaf areas from meshes (parallel)                  │
└─────────────────────────────────────────────────────────────────┘
                       ↓
┌─────────────────────────────────────────────────────────────────┐
│ Output: 2024-02-23_angola_P02_treeLevel_attr_rct_*.csv          │
│ ───────────────────────────────────────────────────────────────  │
│ Contains: plot_id, tile, tree_id, x, y, z, height, DBH,         │
│          crown_radius, in_plot, leaf_area_m2, num_leaves        │
└─────────────────────────────────────────────────────────────────┘

```

### Step 2: MIMICS Parameter Generation

```
┌─────────────────────────────────────────────────────────────────┐
│ Input Files                                                      │
│ ───────────────────────────────────────────────────────────────  │
│ • Tree attributes CSV (from Step 1)                             │
│ • angola_P02_matched_stems.csv (census ↔ RCT mapping)           │
│ • wood_density.csv (species → density lookup)                   │
└─────────────────────────────────────────────────────────────────┘
                       ↓
┌─────────────────────────────────────────────────────────────────┐
│ For Tree 282 in Tile_0_2 (create_parameters.py)                │
│ ───────────────────────────────────────────────────────────────  │
│ 1. Load: angola_p02_raycloud_0_2_trees_282_info.txt             │
│ 2. Parse segments: branch_order, x, y, z, diameter, length      │
│ 3. Calculate branch angles from segment endpoints               │
│ 4. Fit MIMICS PDFs (18 types) using KL divergence              │
│ 5. Get tree attributes: height=12.5m, DBH=27.6cm               │
│ 6. Get leaf area: 85.3 m²                                      │
│ 7. Get wood density: "Burkea africana" → 0.68 g/cm³            │
└─────────────────────────────────────────────────────────────────┘
                       ↓
┌─────────────────────────────────────────────────────────────────┐
│ Parameter Combinations Generation                                │
│ ───────────────────────────────────────────────────────────────  │
│ For frequency in [0.43, 1.2] GHz:                               │
│   For canopy_density in [0.015, 0.06, 0.24, 0.48, 0.72]:       │
│     For soil_moisture in [0.25, 0.5, 0.75]:                     │
│       Create complete parameter set with:                       │
│       • Sensor: frequency, angle                                │
│       • Ground: soil_moisture, roughness, soil composition      │
│       • Tree: trunk dims, crown thickness, canopy density       │
│       • Branches 1-4: length, diameter, density, PDF type      │
│       • Leaves: density calculated from leaf area               │
│                                                                  │
│ → 30 parameter combinations per tree                            │
└─────────────────────────────────────────────────────────────────┘
                       ↓
┌─────────────────────────────────────────────────────────────────┐
│ Output: mimics_inputs_angola_P02.csv                            │
│ ───────────────────────────────────────────────────────────────  │
│ Each row = complete MIMICS parameter set for one scenario       │
│ Columns: tree metadata + 40+ MIMICS model parameters           │
└─────────────────────────────────────────────────────────────────┘

```
### Step 3: MIMICS Model Execution

```
┌─────────────────────────────────────────────────────────────────┐
│ Input: mimics_inputs_angola_P02.csv (run_model.py)              │
└─────────────────────────────────────────────────────────────────┘
                       ↓
┌─────────────────────────────────────────────────────────────────┐
│ For Each Parameter Set (Parallel Execution)                     │
│ ───────────────────────────────────────────────────────────────  │
│ 1. Create temporary MIMICS workspace                            │
│ 2. Configure input files:                                       │
│    • sensor.input (frequency, angle)                           │
│    • ground.input (soil properties)                            │
│    • trunk_and_gross_canopy.input (tree structure)             │
│    • branch_primary.input through branch_5th.input             │
│    • leaf.input (leaf density)                                 │
│ 3. Execute MIMICS executable                                    │
│ 4. Parse output: extract σ₀ backscatter values                 │
│ 5. Clean up temporary files                                     │
└─────────────────────────────────────────────────────────────────┘
                       ↓
┌─────────────────────────────────────────────────────────────────┐
│ Output: mimics_outputs_angola_P02.csv                           │
│ ───────────────────────────────────────────────────────────────  │
│ Contains: input parameters + MIMICS results (σ₀ values)         │
│ Ready for analysis, validation, and radar simulation            │
└─────────────────────────────────────────────────────────────────┘

═══════════════════════════════════════════════════════════════════
```

## PERFORMANCE COMPARISON

```
┌───────────────────────┬──────────────────┬──────────────────────┐
│ Metric                │ Version 1.0       │ Version 2.0         │
├───────────────────────┼──────────────────┼──────────────────────┤
│ Memory Usage          │ High (~5GB)      │ Low (~500MB)         │
│ Disk Usage            │ +2GB redundant   │ +Intermediate CSVs   │
│ Processing Time       │ 10-15 min        │ 3-8 min (parallel)   │
│ Step 1 Focus          │ Segment export   │ Tree attributes ⭐   │
│ Step 2 Processing     │ Sequential       │ Parallel ⚡          │
│ Step 3 Features       │ Basic parallel   │ Enhanced ⚡          │
│ Spatial Filtering     │ ❌ No            │ ✅ Yes (Step 1)      │
│ Leaf Area Calc        │ ❌ No            │ ✅ Yes (Step 1)      │
│ On-demand Loading     │ ❌ No            │ ✅ Yes (Step 2)      │
│ Retry Logic           │ ❌ Limited       │ ✅ Robust (Step 3)   │
│ Progress Tracking     │ ❌ Basic         │ ✅ Detailed          │
│ Error Recovery        │ ❌ Difficult     │ ✅ Easy              │
│ Checkpoints           │ ✅ After steps   │ ✅ After steps       │
│ Reusability           │ ⚠️ Limited       │ ✅ High (modular)    │
│ Debugging             │ ⚠️ Difficult     │ ✅ Easy              │
│ Scalability           │ ⚠️ Limited       │ ✅ Excellent         │
│ Code Maintenance      │ ⚠️ 3 scripts     │ ✅ 3 modular scripts │
└───────────────────────┴──────────────────┴──────────────────────┘

```
### Key Improvements Summary:
- **Step 1**: New comprehensive attribute extraction (not just segment export)
- **Step 2**: Corrected structural attributes definition. Refactored for parallel processing and memory efficiency
- **Step 3**: Enhanced with retry logic, better error handling, progress tracking

═══════════════════════════════════════════════════════════════════


### Dependencies Between Steps

```
RCT Outputs ──► Step 1 ──► Tree Attributes CSV
                             │
Census Data ─────────────────┼──► Step 2 ──► MIMICS Input CSV
Wood Density ────────────────┘                 │                
                                               │
MIMICS Executable ─────────────────────────────┼──► Step 3 ──► Results CSV
MIMICS Templates ──────────────────────────────┘
```

### Error Recovery and Debugging

- **Step 1 fails**: Check RCT file paths, scan positions, leaf meshes
- **Step 2 fails**: Verify matched stems, check individual tree info files
- **Step 3 fails**: Validate MIMICS installation, check parameter ranges
- **Partial failures**: Each step preserves progress, can restart from checkpoints

## FILE NAMING PATTERNS

### Step 1: RCT Extraction Files (Input)
```
RCT Raw Outputs:
  {country}_{plot}_raycloud_{tile_x}_{tile_y}_trees_info.txt
  {country}_{plot}_raycloud_{tile_x}_{tile_y}_trees.txt
  {country}_{plot}_raycloud_{tile_x}_{tile_y}_leaves.ply
                ↓
  angola_p02_raycloud_0_2_trees_info.txt
  angola_p02_raycloud_0_2_trees.txt
  angola_p02_raycloud_0_2_leaves.ply

Individual Tree Files:
  {country}_{plot}_raycloud_{tile_x}_{tile_y}_treesplit/
    {country}_{plot}_raycloud_{tile_x}_{tile_y}_trees_{tree_id}_info.txt
                ↓
  angola_p02_raycloud_0_2_treesplit/
    angola_p02_raycloud_0_2_trees_282_info.txt
```

### Step 1: Tree Attributes Output
```
{date}_{country}_{plot}_treeLevel_attr_rct_tile{len}m_overlap{overlap}m.csv
                ↓
2024-02-23_angola_P02_treeLevel_attr_rct_tile25m_overlap5m.csv
```

### Step 2: Additional Input Files
```
Matched Stems (census ↔ RCT mapping):
  {country}_{PLOT}_matched_stems.csv
                ↓
  angola_P02_matched_stems.csv

Wood Density Lookup:
  wood_density.csv
```

### Step 2: MIMICS Parameters Output
```
mimics_inputs_{date}_{country}_{plot}.csv  (if date available)
mimics_inputs_{country}_{plot}.csv         (if no date)
                ↓
mimics_inputs_2024-02-23_angola_P02.csv
mimics_inputs_angola_P02.csv
```

### Step 3: MIMICS Results Output
```
mimics_outputs_{country}_{plot}.csv
                ↓
mimics_outputs_angola_P02.csv

Individual Results (optional):
  {country}_{date}_{plot}_{tile_coords}_Tree_{tree_id}/
                ↓
  angola_2024-02-23_P02_0_2_Tree_282/
```
