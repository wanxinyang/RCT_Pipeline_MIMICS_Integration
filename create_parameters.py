import pandas as pd
import glob
import numpy as np
from scipy import stats
import os

def mimics_pdf(theta_deg, pdf_type):
    """Calculate MIMICS PDFs"""
    theta_rad = np.deg2rad(theta_deg)
    
    if pdf_type == 1:  # Uniform
        return np.ones_like(theta_deg) / 180.0
    elif pdf_type == 2:  # (sin(2θ))^4
        return (np.sin(2 * theta_rad))**4
    elif pdf_type == 3:  # (sin(θ))^2
        return (np.sin(theta_rad))**2
    elif pdf_type == 4:  # (sin(θ))^3
        return (np.sin(theta_rad))**3
    elif pdf_type == 5:  # (sin(θ))^4
        return (np.sin(theta_rad))**4
    elif pdf_type == 6:  # (sin(θ))^5
        return (np.sin(theta_rad))**5
    elif pdf_type == 7:  # (sin(θ))^6
        return (np.sin(theta_rad))**6
    elif pdf_type == 8:  # (sin(θ))^7
        return (np.sin(theta_rad))**7
    elif pdf_type == 9:  # (sin(θ))^8
        return (np.sin(theta_rad))**8
    elif pdf_type == 10:  # (sin(θ))^9
        return (np.sin(theta_rad))**9
    elif pdf_type == 11:  # (sin(θ-30°))^9
        return (np.sin(theta_rad - np.deg2rad(30)))**9
    elif pdf_type == 12:  # (sin(θ+60°))^9
        return (np.sin(theta_rad + np.deg2rad(60)))**9
    elif pdf_type == 13:  # (sin(θ))^20
        return (np.sin(theta_rad))**20
    elif pdf_type == 14:  # (sin(θ+30°))^20
        return (np.sin(theta_rad + np.deg2rad(30)))**20
    elif pdf_type == 17:  # (sin(θ+64°))^16
        return (np.sin(theta_rad + np.deg2rad(64)))**16
    elif pdf_type == 18:  # (sin(θ+41°))^8
        return (np.sin(theta_rad + np.deg2rad(41)))**8
    else:
        return np.zeros_like(theta_deg)

def calculate_branch_angles(df):
    """Calculate branch direction angles from first to last segment"""
    df = df.copy()
   
    # initialize angle columns
    df['branch_angle_from_vertical'] = np.nan
   
    # group by each complete branch
    for (tree_id, branch_id, branch_order), group in df.groupby(['tree_id', 'branch', 'branch_order']):
        # Sort by position in branch
        group = group.sort_values('pos_in_branch')

        #  ensures that a branch has at least 2 segments
        if len(group) >= 2:
            # Get first and last segments
            start_segment = group.iloc[0]
            end_segment = group.iloc[-1]
           
            # calculate overall branch direction vector
            dx = end_segment['x'] - start_segment['x']
            dy = end_segment['y'] - start_segment['y']
            dz = end_segment['z'] - start_segment['z']
           
            # Calculate length (3D Euclidean distance formula - calculates the straight-line distance between two points in 3D space)
            length = np.sqrt(dx**2 + dy**2 + dz**2)
           
            if length > 0:    
                # Branch angle from vertical
                vertical_angle = np.arccos(abs(dz) / length) * 180 / np.pi
               
                # Apply angles to all segments in this branch
                for idx in group.index:
                    df.loc[idx, 'branch_angle_from_vertical'] = vertical_angle
   
    return df

def calculate_pdf_matches(df):
    """Calculate best MIMICS PDF match for each species/branch_order combination"""
    
    mimics_types = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 17, 18]
    theta = np.linspace(0, 180, 1000)
    
    pdf_results = {}
    
    # group.. 
    for (census_species, branch_order), group in df.groupby(['census_species', 'branch_order']):

        #  ensure that there are at least 10 branch angle measurements before trying to fit a PDF
        if len(group) < 10 or branch_order == 0 or pd.isna(census_species):
            continue
        
        # Get branch-level data (one angle per branch)
        branch_data = group.groupby('branch').agg({
            'branch_angle_from_vertical': 'first',
            'volume': 'sum'
        }).reset_index()
        
        angles = branch_data['branch_angle_from_vertical'].values
        weights = branch_data['volume'].values
        
        # Filter valid data
        valid_mask = (np.isfinite(angles) & np.isfinite(weights) & (weights > 0))
        angles = angles[valid_mask]
        weights = weights[valid_mask]
        
        if len(angles) < 5:
            continue
        
        # Create weighted KDE
        weights = weights / np.sum(weights)
        weight_multiplier = 1000
        weight_counts = np.round(weights * weight_multiplier).astype(int)
        weight_counts = np.maximum(weight_counts, 1)
        
        weighted_angles = []
        for angle, count in zip(angles, weight_counts):
            weighted_angles.extend([angle] * count)
        
        if len(weighted_angles) > 1:
            kde = stats.gaussian_kde(weighted_angles)
            data_density = kde(theta)
            data_density = data_density / np.trapezoid(data_density, dx=theta[1]-theta[0])
        else:
            continue
        
        best_kl = np.inf
        best_pdf = None
        
        # then test each PDF type
        for pdf_type in mimics_types:
            mimics_pdf_values = mimics_pdf(theta, pdf_type)
            mimics_pdf_values = np.maximum(mimics_pdf_values, 0)
            integral = np.trapezoid(mimics_pdf_values, dx=theta[1]-theta[0])
            if integral > 0:
                mimics_pdf_normalized = mimics_pdf_values / integral
            else:
                continue
            
            # KL divergence
            data_safe = np.maximum(data_density, 1e-10)
            mimics_safe = np.maximum(mimics_pdf_normalized, 1e-10)
            kl_div = np.trapezoid(data_safe * np.log(data_safe / mimics_safe), dx=theta[1]-theta[0])
            
            if np.isfinite(kl_div) and kl_div < best_kl:
                best_kl = kl_div
                best_pdf = pdf_type
        
        # Store results
        pdf_results[(census_species, branch_order)] = {'best_pdf': best_pdf, 'kl_divergence': best_kl}
    
    return pdf_results


def calculate(data_dir, census_data, wood_density_file, min_branch_order=4, max_branch_order=4):
    """
    Process tree data files to extract tree characteristics and merge with census/density data.
    
    Args:
        data_dir: Directory containing CSV files with tree data
        census_data: Path to census data CSV file
        wood_density_file: Path to wood density data CSV file
        min_branch_order: Minimum branch order required for tree inclusion (default: 3)
        max_branch_order: Maximum branch order to process (default: 3)
    
    Returns:
        DataFrame with processed tree data including calculated metrics and merged census info
    """
    
    csv_files = glob.glob(os.path.join(data_dir, "angola_p*_combined.csv"))
    all_trees = []
   
    preserve_columns = ['tree_id', 'plot_id', 'tile_coords', 'file_name', 
                       'height', 'crown_radius', 'DBH', 'leaf_area_m2', 
                       'num_leaves', 'cstart']
   
    # Load census data and wood density data
    census_df = pd.read_csv(census_data)
    wood_density_df = pd.read_csv(wood_density_file)

    
   # loop through each CSV file found by the glob pattern
    for file_path in csv_files:
        df = pd.read_csv(file_path, low_memory=False)

        # Calculate branch angles
        df = calculate_branch_angles(df)
        
        # Calculate geometric properties for each segment
        # Volume: cylinder volume = π * radius² * length
        df['volume'] = 3.14159 * (df['diameter'] / 2) ** 2 * df['segment_length']
        # Surface area: cylinder surface area = 2 * π * radius * length
        df['surface_area'] = 2 * 3.14159 * (df['diameter'] / 2) * df['segment_length']
        
        # Calculate PDF matches for this file's data
        pdf_results = calculate_pdf_matches(df)
        
        # For each tree (identified by filename)
        unique_filenames = df['file_name'].unique()
        for file_name in unique_filenames:
            # Filter data for current tree/file
            tree_data_full = df[df['file_name'] == file_name]
            tree_row = tree_data_full.iloc[0]
            
            # Check if tree has required minimum branch order
            max_order_in_tree = tree_data_full['branch_order'].max()
            if max_order_in_tree < min_branch_order:
                continue  # Skip this tree
            
            tree_data = {'file_name': file_name}
            
            # Copy preserved columns from original data
            for col in preserve_columns:
                if col in df.columns:
                    if col == 'DBH':
                        # Convert DBH from meters to cm
                        tree_data[col] = tree_row[col] * 100
                    else:
                        tree_data[col] = tree_row[col]

            
            # Wood density logic
            # get species 
            tree_data['census_species'] = tree_row.get('census_species', None)
            # Match wood density
            if tree_data['census_species'] and pd.notna(tree_data['census_species']):
                density_match = wood_density_df[wood_density_df['census_species'] == tree_data['census_species']]
                
                if not density_match.empty:
                    tree_data['wood_density_mean'] = density_match.iloc[0]['wood_density_mean']
                else:
                    tree_data['wood_density_mean'] = None
            else:
                tree_data['wood_density_mean'] = None

            
            # Calculate crown height and volume
            crown_height = tree_row['height'] - tree_row['cstart']
            crown_volume = (4/3) * 3.14159 * (tree_row['crown_radius'] ** 2) * crown_height
            tree_data['crown_height'] = crown_height
            tree_data['crown_volume'] = crown_volume

            # Canopy Density
            tree_data['canopy_density'] = 0.015
            # Frequency column
            tree_data['frequency'] = 0.5
            tree_data['angle'] = 30
            tree_data['soil_moisture'] = 0.5
            tree_data['rms_height'] = 1.5
            tree_data['correlation_length'] = 17.5
            tree_data['percent_sand'] = 40
            tree_data['percent_clay'] = 10
            # wood moisture columns
            tree_data['trunk_moisture'] = 0.5
            tree_data['branch_1_moisture'] = 0.5
            tree_data['branch_2_moisture'] = 0.5
            tree_data['branch_3_moisture'] = 0.5
            tree_data['branch_4_moisture'] = 0.5
            tree_data['branch_5_moisture'] = 0.5         
            # wood density columns - get the actual density value
            wood_density_value = tree_data.get('wood_density_mean', 0.5)  # Default to 0.5 if no density found
            tree_data['trunk_dry_density'] = wood_density_value
            tree_data['branch_1_dry_density'] = wood_density_value
            tree_data['branch_2_dry_density'] = wood_density_value
            tree_data['branch_3_dry_density'] = wood_density_value
            tree_data['branch_4_dry_density'] = wood_density_value
            tree_data['branch_5_dry_density'] = wood_density_value

            # Calculate branch order statistics (up to max_branch_order)
            for order in range(0, max_branch_order + 1):
                order_data = tree_data_full[tree_data_full['branch_order'] == order]
                
                if len(order_data) > 0:  # If this order exists for this tree
                    # Total volume for this order
                    tree_data[f'branch_{order}_volume'] = order_data['volume'].sum()
                    
                    # Total surface area for this order
                    tree_data[f'branch_{order}_surface_area'] = order_data['surface_area'].sum()
                    
                    # Branch-level calculations
                    branch_stats_list = []
                    
                    for branch_id, branch_group in order_data.groupby('branch'):
                        # Total length per branch (sum of all segments)
                        branch_length = branch_group['segment_length'].sum()
                        
                        # Volume-weighted diameter per branch
                        branch_volume = branch_group['volume'].sum()
                        if branch_volume > 0:
                            weighted_diameter = np.average(branch_group['diameter'], weights=branch_group['volume'])
                        else:
                            weighted_diameter = branch_group['diameter'].mean()
                        
                        branch_stats_list.append({
                            'branch': branch_id,
                            'length': branch_length,
                            'diameter': weighted_diameter,
                            'volume': branch_volume
                        })
                    
                    branch_stats = pd.DataFrame(branch_stats_list)
                    
                    # Volume-weighted average branch length for this order
                    if branch_stats['volume'].sum() > 0:
                        tree_data[f'branch_{order}_length'] = np.average(branch_stats['length'], weights=branch_stats['volume'])
                    else:
                        tree_data[f'branch_{order}_length'] = branch_stats['length'].mean()
                    
                    # Volume-weighted average branch diameter for this order (convert to cm)
                    if branch_stats['volume'].sum() > 0:
                        tree_data[f'branch_{order}_diameter'] = np.average(branch_stats['diameter'], weights=branch_stats['volume']) * 100
                    else:
                        tree_data[f'branch_{order}_diameter'] = branch_stats['diameter'].mean() * 100
                    
                    # Skip density for branch order 0 (this is the trunk)
                    if order != 0:
                        # Number density (number of branches / crown volume)
                        unique_branches = order_data['branch'].nunique()
                        tree_data[f'branch_{order}_density'] = unique_branches / crown_volume if crown_volume > 0 else 0
                    
                    # Skip PDF for branch order 0
                    if order != 0:
                        # Add PDF results if available
                        species_name = tree_data.get('census_species')
                        if species_name and (species_name, order) in pdf_results:
                            pdf_info = pdf_results[(species_name, order)]
                            tree_data[f'branch_{order}_pdf'] = pdf_info['best_pdf']
                        else:
                            tree_data[f'branch_{order}_pdf'] = None
                        
                else:  # If this order doesn't exist for this tree, set to None
                    tree_data[f'branch_{order}_volume'] = None
                    tree_data[f'branch_{order}_surface_area'] = None
                    tree_data[f'branch_{order}_length'] = None
                    tree_data[f'branch_{order}_diameter'] = None
                    
                    # Skip density and PDF for branch order 0
                    if order != 0:
                        tree_data[f'branch_{order}_density'] = None
                        tree_data[f'branch_{order}_pdf'] = None
            
            all_trees.append(tree_data)
   
    # Add calculated ratios and totals
    for tree in all_trees:
        # Calculate totals
        total_tree_volume = sum(tree.get(f'branch_{order}_volume', 0) or 0 for order in range(0, max_branch_order + 1))
        total_tree_surface_area = sum(tree.get(f'branch_{order}_surface_area', 0) or 0 for order in range(0, max_branch_order + 1))
        
        # Canopy totals (exclude order 0 which is trunk)
        total_canopy_volume = sum(tree.get(f'branch_{order}_volume', 0) or 0 for order in range(1, max_branch_order + 1))
        total_canopy_surface_area = sum(tree.get(f'branch_{order}_surface_area', 0) or 0 for order in range(1, max_branch_order + 1))
        
        trunk_volume = tree.get('branch_0_volume', 0) or 0
        
        # Add totals to tree data
        tree['total_tree_volume'] = total_tree_volume
        tree['total_tree_surface_area'] = total_tree_surface_area
        tree['total_canopy_volume'] = total_canopy_volume
        tree['total_canopy_surface_area'] = total_canopy_surface_area
        
        # Calculate ratios (avoid division by zero)
        if total_tree_volume > 0:
            tree['tree_sa_to_volume_ratio'] = total_tree_surface_area / total_tree_volume
            tree['trunk_to_tree_volume_ratio'] = trunk_volume / total_tree_volume
            tree['canopy_to_tree_volume_ratio'] = total_canopy_volume / total_tree_volume
            
            # Order volume to tree volume ratios
            for order in range(0, max_branch_order + 1):
                order_volume = tree.get(f'branch_{order}_volume', 0) or 0
                tree[f'order_{order}_to_tree_volume_ratio'] = order_volume / total_tree_volume
        else:
            tree['tree_sa_to_volume_ratio'] = None
            tree['trunk_to_tree_volume_ratio'] = None
            tree['canopy_to_tree_volume_ratio'] = None
            for order in range(0, max_branch_order + 1):
                tree[f'order_{order}_to_tree_volume_ratio'] = None
        
        # Order surface area to order volume ratios
        for order in range(0, max_branch_order + 1):
            order_volume = tree.get(f'branch_{order}_volume', 0) or 0
            order_surface_area = tree.get(f'branch_{order}_surface_area', 0) or 0
            if order_volume > 0:
                tree[f'order_{order}_sa_to_volume_ratio'] = order_surface_area / order_volume
            else:
                tree[f'order_{order}_sa_to_volume_ratio'] = None
   
    # Filter out trees with missing PDF data for any branch order
    filtered_trees = []
    for tree in all_trees:
        has_all_pdfs = True
        for order in range(1, max_branch_order + 1):
            if tree.get(f'branch_{order}_pdf') is None:
                has_all_pdfs = False
                break
        
        if has_all_pdfs:
            filtered_trees.append(tree)
   
    # Calculate crown height and crown volume
    df_result = pd.DataFrame(filtered_trees)


    # Create multiple datasets with different parameter values
    frequency_values = [0.43]
    soil_moisture_values = [0.5]  

    expanded_data = []
    for freq_val in frequency_values:
        for soil_moist_val in soil_moisture_values:
            df_copy = df_result.copy()
            df_copy['frequency'] = freq_val
            df_copy['soil_moisture'] = soil_moist_val
            expanded_data.append(df_copy)
    
    # Combine all datasets
    df_result = pd.concat(expanded_data, ignore_index=True)
    
    # Rename columns - adjust mappings here as needed
    column_rename_map = {
        'tree_id': 'tree_id',
        'plot_id': 'plot_id',
        'tile_coords': 'tile_coords',
        'file_name': 'file_name',
        'height': 'height',
        'crown_radius': 'crown_radius',
        'branch_0_volume': 'trunk_volume',       
        'branch_0_surface_area': 'trunk_surface_area',       
        'branch_0_length': 'trunk_length',
        'DBH': 'trunk_diameter',
        'leaf_area_m2': 'leaf_density',
        'num_leaves': 'num_leaves',
        'cstart': 'cstart',
        'census_species': 'census_species',
        'wood_density_mean': 'wood_density',
        'crown_height': 'crown_thickness',
        'crown_volume': 'crown_volume',
        'frequency': 'frequency',
        'angle': 'angle',
        'soil_moisture': 'soil_moisture',
        'correlation_length': 'correlation_length',
        'percent_sand': 'percent_sand',
        'percent_clay': 'percent_clay',
    }
    
    df_result = df_result.rename(columns=column_rename_map)
   
    return df_result


if __name__ == "__main__":
    tree_data = calculate(
        data_dir=r"/home/ucfargt@ad.ucl.ac.uk/Documents/mimics/segment_data",
        census_data=r"/home/ucfargt@ad.ucl.ac.uk/Documents/mimics/census_data.csv",
        wood_density_file=r"/home/ucfargt@ad.ucl.ac.uk/Documents/mimics/wood_density.csv",
        min_branch_order=4, 
        max_branch_order=4
    )
    
    tree_data.to_csv(r"/home/ucfargt@ad.ucl.ac.uk/Documents/mimics/model_input.csv", index=False)
    print(f"Found {len(tree_data)} trees")
