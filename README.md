# To-do
-  Add input file safety logic to protect when cancelling a run
-  add logic that adjusts mimics config file for which branch order enabled/disabled so you dont have to do it manually

# RCT-Pipeline to MIMICS Integration
Three-script workflow for processing 3D tree structure data from RCT-pipeline (https://github.com/wanxinyang/rct-pipeline
) and running radiative scattering simulations using the Michigan Microwave Canopy Scattering Model (MIMICS; https://codeocean.com/capsule/8976769/tree/v1).

## Purpose
i. Extract segment data from Quantitative Structure Models (QSMs) into a csv (per plot and tile)  
ii. Calculate tree parameters suitable for MIMICS  
iii. Input TLS derived parameters into the MIMICS, run the model, and output a backscatter reponse dataset.  

### Recommended Directory Structure

            mimics/  
            ├── model/  
            │   ├── code
            │   ├── data
            │   └── ...  
            └── get_tree_data.py # script 1
            └── create_parameters.py # script 2
            └── run_model.py # script 3
            └── census_data.csv
            └── wood_density.csv
            └── segment_data # script 1 output
            └── model_input.csv # script 2 output
            └── model_output.csv # script 3 output

## Script 1: get_tree_data.py
### Overview
The first script to run is get_tree_data.py, which creates CSV files containing segment data extracted from the Quantitative Structure Models (QSMs). The output files are grouped per plot and tile.

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
3. Census file containing addtional quality assessment (see census_data.csv) and the suitability coloumn (trees marked as 'good') Plot IDs, tile coordinates, and tree IDs Addtionally, the height at which the crown starts (c_start) should also be recorded for later calculations. 

You will have to update this section of the script accordingly:  

            if __name__ == "__main__":
                tree_files = "/path/to/rct_extraction"  
                census_data = "/path/to/census_data.csv" 

### Output
Creates angola_{plot_id}_{tile_coords}_combined.csv files in ./segment_data/ directory in the same location where the Python script is saved and run from.

## Script 2: create_parameters.py
Processes QSM segment data to create MIMICS model input parameters. Calculates branch angles, matches optimal PDF distributions for species/branch combinations using KL divergence, and generates volume/surface area/density statistics per branch order.

### Prerequisites
Required Files:
1. Output CSV files from Script 1 (angola_*_combined.csv)  
2. census_data.csv - identification data  
3. wood_density.csv - wood density values by species  


### Configuration 
This script is the primary config file and produces the input parameter files.  

First assign the files accordingly.

            if __name__ == "__main__":
                tree_data = calculate(
                    data_dir=r"/path/to/segment_data",
                    census_data=r"/path/to/census_data.csv",
                    wood_density_file=r"/path/to/wood_density.csv",
                    min_branch_order=4, 
                    max_branch_order=4
                )
                
Min and max branch order must also be defined. 
1. min_branch_order: Minimum branch order required for trees to be included in the simulation (default: 4).  
2. max_branch_order: Maximum branch order for calculating statistics (default: 4).

Note: code currently lacks dynamic logic to automatically adjust the mimics.configuration.input file based on the branch orders of individual tree files. This creates a problem because different trees can have different branch order distributions. To work around this, the configuration file is manually set to simulate specific branch orders (for example, orders 1 through 4), and then we filter the input dataset to only include trees that have at least the minimum branch order (in this case, 4). This ensures that every tree in the filtered dataset has sufficient structural complexity to match what the configuration file expects.  

### Data-derived 
The values derived from the tree files are: 
1. Branch angle probability density functions (PDFs) calculated at the species-level, not tree-level. All trees of the same species share the same PDF values. assumption is that branching architecture is a species characteristic.  
3. Branch order diameter, length, volume, surface area, density  
4. Trunk, length, diameter, surface area and volume
5. Canopy volume and height
7. Volume ratios between branch orders  

### Hard-coded
The following parameters cannot be derived from RCT data and are set here:

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


### creating multiple datasets to sweep. 
To create a large dataset that sweeps across various parameter combinations, edit the parameter value lists in this section of the code. the loops will generate all possible combinations of the specified parameters (e.g., frequency, soil moisture, and any additional input parameters). Creating dataset copies for each combination.

    # Create multiple datasets with different parameter values
    frequency_values = [0.43, 0.5, 1.2, 5.1]
    soil_moisture_values = [0.5]  
    # input_parameter_values = [1, 2, 3, 4, 5] **(example)**

    expanded_data = []
    for freq_val in frequency_values:
        for soil_moist_val in soil_moisture_values:
            for input_parameter_val in input_parameter_values: **(example)**      
                df_copy = df_result.copy()
                df_copy['frequency'] = freq_val
                df_copy['soil_moisture'] = soil_moist_val
                # df_copy['input_parameter'] = input_parameter_val **(example)**
                expanded_data.append(df_copy)

### Output
model_input_data.csv 

## Script 3. run_model.py
### 5.1. Overview 
Takes a CSV row (from model_input_data.csv) and writes the parameter values into the correct MIMICS model input files at the right line numbers and formatting. This uses multi_threading. 

### Prerequisites  

model_input_data.csv from Script 2  
MIMICS model files as described earlier  

### Output  
model_output.csv

## Analysis

Depending on what was chosen to sweep in the creating multiple datasets section, The sav_analysis script allows to explore these results. Currently this includes frequency, soil moisture and canopy density.  In addititon, the x_param  e.g., tree_sa_to_volume_ratio, tree_total_volume, etc. can be changed and we can also choose the backscatter mechanism  e.g., backscatter_type='total', ground_trunk, direct_ground etc. 

                        x_param='tree_sa_to_volume_ratio', backscatter_type='total', frequency=0.43, angle=30, soil_moisture=0.5, canopy_density=0.015)












