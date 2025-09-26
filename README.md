# RCT-Pipeline Processing Scripts

This is a set of four scripts designed to follow the rct-pipeline workflow. The purpose is to:  

i. Output a segment level csv dataset for each tile and plot.  
ii. Create an input parameter dataset suitable for the mimics model  
iii. Inputs the parameters into the mmimics model  
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

3. Reference file: A quality assessment file (see matched_stems_picked.csv) containing:
Tree suitability ratings (must include trees marked as 'good')
Plot IDs, tile coordinates, and tree IDs
Addtionally, the height at which the crown starts should also be recorded for later calculations. 

Update these paths in script:
tree_files_dir = "/path/to/rct_extraction"
matched_stems_file = "/path/to/matched_stems_picked.csv"

### Output
Creates angola_{plot_id}_{tile_coords}_combined.csv files in ./output/ directory. The "output" directory in the same location where the Python script is saved and run from.


## 4. create_parameters.py
This next script transforms raw segment-level tree data (Step 1 Output) into tree-level input parameter sets suitable for the mimics model. 

### 4.1. Running
The min_branch_order and max_branch_order parameters control which trees are included in the analysis and what branch data is calculated for each tree. This is for both quality control - ensures all trees have suffienent branching - and consistency - All trees in the final dataset have the same range of branch orders.

min_branch_order (Quality Filter) sets the minimum branch order requirement for trees to be included in the dataset. By default min_branch_order=3, meaning, only trees that have at least branch order 3 will be included. 
max_branch_order (Data Scope) determines how many branch orders to calculate statistics for. 

#### 4.2 Parameters 

Some parameters cannot be derived from the rct data.  These values are simply hard-coded and include the following: 

            tree_data['canopy_density'] = 0.015 # stocking density
            
            tree_data['frequency'] = 0.5
            tree_data['angle'] = 30

            tree_data['soil_moisture'] = 0.5
            tree_data['rms_height'] = 1.5
            tree_data['correlation_length'] = 17.5
            tree_data['percent_sand'] = 40
            tree_data['percent_clay'] = 10
            
            tree_data['trunk_moisture'] = 0.5
            tree_data['branch_1_moisture'] = 0.5
            tree_data['branch_2_moisture'] = 0.5
            tree_data['branch_3_moisture'] = 0.5
            tree_data['branch_4_moisture'] = 0.5
            tree_data['branch_5_moisture'] = 0.5
            
            wood_density_value = tree_data.get('wood_density_mean', 0.5) 
            tree_data['trunk_dry_density'] = wood_density_value
            tree_data['branch_1_dry_density'] = wood_density_value
            tree_data['branch_2_dry_density'] = wood_density_value
            tree_data['branch_3_dry_density'] = wood_density_value
            tree_data['branch_4_dry_density'] = wood_density_value
            tree_data['branch_5_dry_density'] = wood_density_value

Wood density values are provided through the wood_density.csv 


#### 4.2.2. data derived parameters

The following parameters are calculated or directly pulled using data from the rct-pipeline files

Branch Angle Probability Density Functions
Branch Diameter
Branch Length 

#### 4.2.3. Additional Calculations

Trunk Surface Area 
Trunk Volume 


## 5. set_parameters.py

### 5.1. Overview 
This script reads parameter values from a CSV file and updates the corresponding input files with proper formatting.

## 6. run_model.py











