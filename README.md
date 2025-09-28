# RCT-Pipeline Processing Scripts

This is a set of four scripts designed to follow the rct-pipeline workflow. 

The purpose is to:  

i. Output a segment level csv dataset for each tile and plot.  
ii. Create an input parameter dataset suitable for the mimics model  
iii. Input the parameters into the mmimics model  
iv. Run mimics and output a results dataset including all input parameters and radar backscatter reponse.  


## Script 1: get_tree_data.py
### Overview
The first script to run is get_tree_data.py, which creates CSV files containing segment data extracted from Quantitative Structure Models (QSMs). 

### Prerequisites

Before running this script, ensure the rct-pipeline has been used to produce:

Required Directory Structure:

            rct_extraction/  
            ├── angola_p02_raycloud_-2_1_treesplit/  
            │   ├── angola_p02_raycloud_-2_1_trees_12_info.txt  
            │   ├── angola_p02_raycloud_-2_1_segmented_12_leaves.ply  
            │   └── ...  
            └── [other plot/tile directories]/  


Required Files:
1. Tree info files: *_trees_*_info.txt files containing structural attributes of each tree QSM

2. Leaf files: *_leaves.ply files containing 3D mesh data for leaf area calculations

3. Census file containing addtional quality assessment (see census_data.csv) and the suitability coloumn (trees marked as 'good')
Plot IDs, tile coordinates, and tree IDs
Addtionally, the height at which the crown starts (c_start) should also be recorded for later calculations. 

Update these paths in script:  
tree_files_dir = "/path/to/rct_extraction"  
matched_stems_file = "/path/to/census_data.csv"  

### Output
Creates angola_{plot_id}_{tile_coords}_combined.csv files in ./output/ directory. The "output" directory in the same location where the Python script is saved and run from.



## Script 2: create_parameters.py
Processes QSM segment data to create MIMICS model input parameters. Calculates branch angles, matches optimal PDF distributions for species/branch combinations using KL divergence, and generates volume/surface area/density statistics per branch order.

### Prerequisites
1. Output CSV files from Script 1 (angola_*_combined.csv)  
2. census_data.csv - identification data  
3. wood_density.csv - wood density values by species  


### Configuration 

min_branch_order: Minimum branch order required for trees to be included (default: 3)  
max_branch_order: Maximum branch order for calculating statistics (default: 3) # note that not all trees may have the same max order so best tto keep these identical for now (fix this later)  

### Hard-coded
The following parameters cannot be derived from RCT data and are set as constants:

            tree_data['canopy_density'] = 0.015 # stocking density
            tree_data['frequency'] = 0.5 # radar frequency  
            tree_data['angle'] = 30 # incidence angle
            tree_data['soil_moisture'] = 0.5
            tree_data['rms_height'] = 1.5
            tree_data['correlation_length'] = 17.5
            tree_data['percent_sand'] = 40
            tree_data['percent_clay'] = 10
            
            # Moisture content (per branch order)
            tree_data['trunk_moisture'] = 0.5
            tree_data['branch_1_moisture'] = 0.5   # ... through branch_5_moisture

            # Wood density (derived from wood_density.csv or default 0.5)

            wood_density_value = tree_data.get('wood_density_mean', 0.5) 
            tree_data['trunk_dry_density'] = wood_density_value
            tree_data['branch_1_dry_density'] = wood_density_value  # ... through branch_5_dry_density

### Data-derived 

Branch angle probability density functions (PDFs)  
Branch diameter, length, volume, surface area  
Trunk surface area and volume  
Volume ratios between branch orders  

### Output
model_input_data.csv containing tree-level parameters with columns for each branch order's volume, surface area, length, diameter, density, and optimal PDF type.
            

## Script 3. set_parameters.py  
### 5.1. Overview 
Takes a CSV row (from model_input_data.csv) and writes the parameter values into the correct MIMICS model input files at the right line numbers and formatting.

### Prerequisites  

model_input_data.csv from Script 2
MIMICS input files in ./model/data/ directory:

### Output  
Updated MIMICS input files ready for model execution with parameters from the specified CSV row.



## Script 4. run_model.py











