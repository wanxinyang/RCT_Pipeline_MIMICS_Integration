
"""
Set MIMICS parameters in input files from CSV
"""
import os
import pandas as pd
import re
from pathlib import Path

# Define paths
base_dir = Path("D:/mimics")
data_dir = base_dir / "model" / "data"

# Parameter definitions: (filename, template_line, value_line)
PARAMETERS = {
    'frequency': ('sensor.input', 13, 14),
    'angle': ('sensor.input', 20, 21),

    'soil_moisture': ('ground.input', 22, 23),
    'rms_height': ('ground.input', 29, 30),
    'correlation_length': ('ground.input', 36, 37),
    'percent_sand': ('ground.input', 45, 46),
    'percent_clay': ('ground.input', 50, 51),
    
    'trunk_moisture': ('trunk_and_gross_canopy.input', 14, 15),
    'trunk_dry_density': ('trunk_and_gross_canopy.input', 21, 22),
    'canopy_density': ('trunk_and_gross_canopy.input', 28, 29),
    'crown_thickness': ('trunk_and_gross_canopy.input', 35, 36),
    'trunk_diameter': ('trunk_and_gross_canopy.input', 42, 43),
    'trunk_length': ('trunk_and_gross_canopy.input', 49, 50),

    'branch_1_moisture': ('branch_primary.input', 14, 15),
    'branch_1_dry_density': ('branch_primary.input', 21, 22),
    'branch_1_length': ('branch_primary.input', 28, 29),
    'branch_1_diameter': ('branch_primary.input', 35, 36),
    'branch_1_density': ('branch_primary.input', 42, 43),
    'branch_1_pdf': ('branch_primary.input', 53, 54),

    'branch_2_moisture': ('branch_secondary.input', 14, 15),
    'branch_2_dry_density': ('branch_secondary.input', 21, 22),       
    'branch_2_length': ('branch_secondary.input', 28, 29),
    'branch_2_diameter': ('branch_secondary.input', 35, 36),
    'branch_2_density': ('branch_secondary.input', 42, 43),
    'branch_2_pdf': ('branch_secondary.input', 53, 54),

    'branch_3_moisture': ('branch_3rd.input', 14, 15),
    'branch_3_dry_density': ('branch_3rd.input', 21, 22),       
    'branch_3_length': ('branch_3rd.input', 28, 29),
    'branch_3_diameter': ('branch_3rd.input', 35, 36),
    'branch_3_density': ('branch_3rd.input', 42, 43),
    'branch_3_pdf': ('branch_3rd.input', 53, 54),

    'branch_4_moisture': ('branch_4th.input', 14, 15),
    'branch_4_dry_density': ('branch_4th.input', 21, 22),       
    'branch_4_length': ('branch_4th.input', 28, 29),
    'branch_4_diameter': ('branch_4th.input', 35, 36),
    'branch_4_density': ('branch_4th.input', 42, 43),
    'branch_4_pdf': ('branch_4th.input', 53, 54),

    'branch_5_moisture': ('branch_5th.input', 14, 15),
    'branch_5_dry_density': ('branch_5th.input', 21, 22),       
    'branch_5_length': ('branch_5th.input', 28, 29),
    'branch_5_diameter': ('branch_5th.input', 35, 36),
    'branch_5_density': ('branch_5th.input', 42, 43),
    'branch_5_pdf': ('branch_5th.input', 53, 54),
    
    'leaf_density': ('leaf.input', 41, 42),

}

def get_exact_formatting_from_template(filename, template_line_num, value_line_num):
    """read MIMICS input files to understand  exact formatting requirements and creates format strings for parameter replacement."""
    
    # create path to mimics input file
    filepath = data_dir / filename

    # opens file and read all lines into a list
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # get the specific lines
    template_line = lines[template_line_num - 1]
    value_line = lines[value_line_num - 1]
    
    # Handle PDF parameters (integer with II pattern)
    # check if integar parameter
    if 'II' in template_line:
        ii_match = re.search(r'II+', template_line)
        if ii_match:
            start_pos = ii_match.start()
            end_pos = ii_match.end()
            field_width = len(ii_match.group(0))
            replacement = f"{{value:{field_width}d}}"
            format_string = (value_line[:start_pos] + replacement + value_line[end_pos:])
            return format_string
    
    # Handle float parameters (XX.xxx pattern)
    # Finds all float patterns like 'XX.xxx', 'XXX.xx', etc.
    template_pattern = r'X+\.\w+'
    template_matches = list(re.finditer(template_pattern, template_line))
    
    # If no patterns found, return None (invalid template)
    if not template_matches:
        return None
    
    # Start with the original value line and track replacements
    result_line = value_line
    replacements_made = 0
    
    # Process patterns (reverse order to avoid shifting issue), skip increment setting if 3 patterns exist 
    for i, match in enumerate(reversed(template_matches)):
        if len(template_matches) == 3 and i == 0:
            continue
        if replacements_made < 2:
            start_pos = match.start()
            end_pos = match.end()
            template_pattern_text = match.group(0)
            
            # Get decimal places from template ('XX.X' = 1 decimal place, 'XXX.xxx' = 3 decimal places)
            if '.' in template_pattern_text:
                decimal_places = len(template_pattern_text.split('.')[1])
            else:
                decimal_places = 0
            
            # Create format specification ('XX.X' (4 chars, 1 decimal) → "{value:4.1f}")
            field_width = len(template_pattern_text)
            format_spec = f"{field_width}.{decimal_places}f"
            replacement = f"{{value:{format_spec}}}"
            
            # replace the pattern with format placeholder and increment counter (e.g.: "frequency  0.5 GHz" becomes "frequency {value:4.1f} GHz"
            result_line = (result_line[:start_pos] + replacement + result_line[end_pos:])
            replacements_made += 1
    
    return result_line

def set_parameters_from_csv(csv_file, row_index=0):
    """Set parameters from CSV row in input files"""
    df = pd.read_csv(csv_file)
    parameters = df.iloc[row_index].to_dict()
    
    # Group parameters by file 
    files_to_modify = {}
    
    # loop through each parameter from the CSV row
    for param_name, value in parameters.items():
        if param_name in PARAMETERS:
            # Get the file and line numbers where this parameter goes (e.g.,'frequency' → ('sensor.input', 13, 14))
            filename, template_line, value_line = PARAMETERS[param_name]
            
            # initialize list for this file if first time seeing it
            if filename not in files_to_modify:
                files_to_modify[filename] = []
            
            # Get exact formatting from template
            format_string = get_exact_formatting_from_template(filename, template_line, value_line)
            
            if format_string:
                files_to_modify[filename].append((value_line, format_string, value))
    
    # Apply Modifications to Each File
    # Process each MIMICS input file that needs changes
    for filename, modifications in files_to_modify.items():

        # read all lines from the current file into memory
        filepath = data_dir / filename
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Apply modifications for this file
        for value_line, format_string, value in modifications:
            array_index = value_line - 1
            if 0 <= array_index < len(lines):
                # Convert to int for PDF parameters
                if 'd}' in format_string:
                    value = int(value)
                lines[array_index] = format_string.format(value=value)
        
        # Write the modified lines back to the file
        with open(filepath, 'w') as f:
            f.writelines(lines)
    
    print(f"Parameters set from row {row_index}")
    return parameters

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python set_parameters.py <csv_file> [row_index]")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    row_index = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    
    parameters = set_parameters_from_csv(csv_file, row_index)
    print("Parameters set successfully")