#!/usr/bin/env python3
import pandas as pd
import subprocess
import re
import shutil
from pathlib import Path
from multiprocessing import Pool, cpu_count

# paths
base_dir = Path("/home/ucfargt@ad.ucl.ac.uk/Documents/mimics/")
data_dir = base_dir / "model" / "data"
code_dir = base_dir / "model" / "code"
results_dir = base_dir / "model" / "results"

# directory to store outputs
preserved_outputs_dir = base_dir / "preserved_outputs"

# parameter locations
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

def get_format(filename, template_line_num, value_line_num, target_data_dir):
    """read MIMICS input files to understand exact formatting requirements and creates format strings for parameter replacement."""

    # create path to mimics input file
    filepath = target_data_dir / filename
    # opens file and read all lines into a list
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # get the specific lines
    template_line = lines[template_line_num - 1]
    value_line = lines[value_line_num - 1]
    
    # Handle PDF parameters (integer with II pattern)
    if 'II' in template_line:
        ii_match = re.search(r'II+', template_line)
        if ii_match:
            start_pos = ii_match.start()
            end_pos = ii_match.end()
            field_width = len(ii_match.group(0))
            return value_line[:start_pos] + f"{{value:{field_width}d}}" + value_line[end_pos:]
    
    # Handle float parameters (XX.xxx pattern)
    # Finds all float patterns like 'XX.xxx', 'XXX.xx', etc.
    template_pattern = r'X+\.\w+'
    template_matches = list(re.finditer(template_pattern, template_line))
    
    # If no patterns found, return None (invalid template)
    if not template_matches:
        return None
    
    # Start with the original value line and track replacements
    result_line = value_line
    replacements = 0
    
    # Process patterns (reverse order to avoid shifting issue), skip increment setting if 3 patterns exist 
    for i, match in enumerate(reversed(template_matches)):
        if len(template_matches) == 3 and i == 0:
            continue
        if replacements < 2:
            start_pos = match.start()
            end_pos = match.end()
            template_text = match.group(0)
            
            # Get decimal places from template ('XX.X' = 1 decimal place, 'XXX.xxx' = 3 decimal places)
            if '.' in template_text:
                decimal_places = len(template_text.split('.')[1])
            else:
                decimal_places = 0
            
            # Create format specification ('XX.X' (4 chars, 1 decimal) → "{value:4.1f}")
            field_width = len(template_text)
            format_spec = f"{field_width}.{decimal_places}f"
            
            # replace the pattern with format placeholder and increment counter (e.g.: "frequency  0.5 GHz" becomes "frequency {value:4.1f} GHz"            
            result_line = result_line[:start_pos] + f"{{value:{format_spec}}}" + result_line[end_pos:]
            replacements += 1
    
    return result_line

def set_parameters(params, target_data_dir):
    """Set parameters in the input files within the specified data directory"""
    # Create empty dictionary to organize which files need changes
    # Structure will be: {filename: [list of modifications]}
    files_to_modify = {}
    
    # Loops through each parameter in the input dictionary 
    for param_name, value in params.items():
        # check if parameter is one we know how to handle (defined in PARAMETERS dictionary)
        if param_name in PARAMETERS:
            # if good, gets file and line numbers for this parameter
            filename, template_line, value_line = PARAMETERS[param_name]
            # If this is the first parameter for this file, create an empty list for it
            if filename not in files_to_modify:
                files_to_modify[filename] = []
            
            # Gets the exact formatting requirements for this parameter. Returns something like "frequency {value:4.1f} GHz"
            format_string = get_format(filename, template_line, value_line, target_data_dir)
            # If formatting worked, add this modification to the file's list (Stores: (line_number, format_string, actual_value))
            if format_string:
                files_to_modify[filename].append((value_line, format_string, value))
    
    # Second Loop - Apply Changes:
    # Process each MIMICS file that needs changes
    for filename, modifications in files_to_modify.items():
        # Read all lines from the current file into memory
        filepath = target_data_dir / filename
        with open(filepath, 'r') as f:
            lines = f.readlines()
        # Apply each modification for this file
        for value_line, format_string, value in modifications:
            # If this is an integer parameter (PDF), convert to integer
            if 'd}' in format_string:
                value = int(value)
            # Replace the specific line with the formatted value (value_line - 1 converts line number to array index)
            lines[value_line - 1] = format_string.format(value=value)
        
        with open(filepath, 'w') as f:
            f.writelines(lines)

def run_model(target_code_dir):
    """Execute the MIMICS model in the specified code directory"""
    # Simply run the model directly on Linux
    process = subprocess.run(['bash', '-c', f'cd {target_code_dir} && echo "go" | ./mimics1.5'],
                           capture_output=True, text=True)
    return process.returncode == 0

def preserve_outputs(run_number, csv_filename, target_results_dir):
    """Preserve outputs for this run in a directory"""
    # Create main preserved outputs directory if it doesn't exist
    preserved_outputs_dir.mkdir(exist_ok=True)
    
    # Create run-specific directory
    run_dir = preserved_outputs_dir / f"{Path(csv_filename).stem}_run_{run_number:04d}"
    run_dir.mkdir(exist_ok=True)
    
    # Copy all output files from results directory
    if target_results_dir.exists():
        for output_file in target_results_dir.glob("*"):
            if output_file.is_file():
                shutil.copy2(output_file, run_dir / output_file.name)

def parse_backscatter(target_results_dir):
    """Parse backscatter results from the specified results directory"""
    like_file = target_results_dir / "forest_sigma_like.out"
    cross_file = target_results_dir / "forest_sigma_cross.out"
    
    # Component column mappings
    components = {
        'total': (1, 2),
        'gnd_crown_gnd': (3, 4),
        'crown_ground': (5, 6),
        'ground_crown': (7, 8),
        'direct_crown': (9, 10),
        'trunk_ground': (11, 12),
        'ground_trunk': (13, 14),
        'direct_ground': (15, 16)
    }
    
    def parse_file(filename, is_cross=False):
        data = {}
        if not filename.exists():
            return data
            
        with open(filename, 'r') as f:
            lines = f.readlines()
        # Find where data starts (2 lines after header with 'Theta' and 'Total Backscatter')
        data_start = next((i + 2 for i, line in enumerate(lines) 
                          if 'Theta' in line and 'Total Backscatter' in line), -1)
        
        # Exit if no data section found
        if data_start == -1:
            return data
        
        # If data is found, process each data line
        for line in lines[data_start:]:
            parts = line.split() # line.split() splits the line into segments at each whitespace.
            try:
                if len(parts) >= 16:
                    theta = float(parts[0])
                    data[theta] = {}
                    # Extract values for each scattering component
                    for comp_name, (col1, col2) in components.items():
                        # Convert to None if value is missing data flag (-900)
                        val1 = float(parts[col1]) if float(parts[col1]) > -900 else None
                        val2 = float(parts[col2]) if float(parts[col2]) > -900 else None
                        # Label columns based on polarization type
                        if is_cross:
                            data[theta][f'{comp_name}_vh'] = val1
                            data[theta][f'{comp_name}_hv'] = val2
                        else:
                            data[theta][f'{comp_name}_vv'] = val1
                            data[theta][f'{comp_name}_hh'] = val2
            except (ValueError, IndexError):
                continue # Skip problem lines
        return data
    # Parse both output files
    like_data = parse_file(like_file, is_cross=False) # VV/HH polarizations
    cross_data = parse_file(cross_file, is_cross=True) # VH/HV polarizations
    
    # Combine results from both files
    results = []
    all_thetas = set(like_data.keys()) | set(cross_data.keys())

    # Create combined data for each angle
    for theta in sorted(all_thetas):
        row = {'theta': theta}
        if theta in like_data:
            row.update(like_data[theta])
        if theta in cross_data:
            row.update(cross_data[theta])
        results.append(row)
    
    return pd.DataFrame(results)

def run_single_model(args):
    """
    Worker function to run model for one parameter set in an isolated directory.
    
    Args:
        args: tuple of (run_number, row_dict, csv_filename, preserve_flag)
    
    Returns:
        DataFrame with results for this run, or None if run failed
    """
    import tempfile
    
    run_number, row_dict, csv_filename, preserve = args
    
    try:
        # Create temporary isolated directory for this worker
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create isolated copies of model directories
            temp_data_dir = temp_path / "data"
            temp_code_dir = temp_path / "code"
            temp_results_dir = temp_path / "results"
            
            # Copy data and code directories
            shutil.copytree(data_dir, temp_data_dir)
            shutil.copytree(code_dir, temp_code_dir)
            temp_results_dir.mkdir()
            
            print(f"Run {run_number}: Starting...")
            
            # Set parameters in isolated directory
            set_parameters(row_dict, temp_data_dir)
            
            # Run model in isolated directory
            if run_model(temp_code_dir):
                # Parse results from isolated directory
                result = parse_backscatter(temp_results_dir)
                
                # Add input parameters to output
                for param, value in row_dict.items():
                    result[param] = value
                result['run'] = run_number - 1  # 0-indexed
                
                # Preserve outputs if requested
                if preserve:
                    preserve_outputs(run_number, csv_filename, temp_results_dir)
                
                print(f"Run {run_number}: Completed")
                return result
            else:
                print(f"Run {run_number}: Failed")
                return None
                
    except Exception as e:
        print(f"Run {run_number}: Error - {str(e)}")
        return None

def main():
    import sys
    
    # Use default CSV file
    csv_file = base_dir / "model_input.csv"
    
    # Check for preserve option
    preserve = "--preserve-outputs" in sys.argv
    
    # Check for number of processes (default to cpu_count)
    num_processes = cpu_count()
    if "--processes" in sys.argv:
        try:
            idx = sys.argv.index("--processes")
            num_processes = int(sys.argv[idx + 1])
        except (IndexError, ValueError):
            print("Invalid --processes argument, using all CPU cores")
    
    # Load and clean CSV data
    df = pd.read_csv(csv_file)
    print(f"CSV loaded: {df.shape} (rows × columns)")
    df = df.dropna(how='all')  # Remove completely empty rows
    print(f"After removing empty rows: {df.shape}")
    
    # Start batch processing
    print(f"Starting batch run with {len(df)} parameter sets")
    print(f"Using {num_processes} parallel processes")
    if preserve:
        print(f"Output files saved to: {preserved_outputs_dir}")
    
    # Prepare arguments for each worker
    work_items = []
    for i, (_, row) in enumerate(df.iterrows()):
        # Convert row to dict for pickling
        row_dict = row.to_dict()
        work_items.append((i + 1, row_dict, csv_file, preserve))
    
    # Run in parallel using multiprocessing
    with Pool(processes=num_processes) as pool:
        results = pool.map(run_single_model, work_items)
    
    # Filter out None results (failed runs) and collect successful ones
    all_results = [r for r in results if r is not None]
    
    # Save combined results
    if all_results:
        output_file = "model_output.csv"
        combined = pd.concat(all_results, ignore_index=True)
        combined.to_csv(output_file, index=False)
        print(f"\n{'='*60}")
        print(f"Processing complete!")
        print(f"Successful runs: {len(all_results)}/{len(df)}")
        print(f"Results saved to: {output_file}")
        print(f"{'='*60}")
    else:
        print("\nNo successful runs to save.")

if __name__ == "__main__":
    main()