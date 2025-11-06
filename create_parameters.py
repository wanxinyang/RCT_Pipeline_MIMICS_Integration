"""
================================================================================
MIMICS Parameter Generation Script
================================================================================
This script generates MIMICS radar backscatter model input parameters from 
files derived from RayCloudTools (RCT) extracted & reconstructed trees.

Main workflow:
1. Load tree attributes and optional census data
2. Validate plot consistency (single plot only)
3. Filter in-plot trees by DBH and height thresholds
4. For each tree:
   - Load RCT segment data from extraction files
   - Calculate branch statistics and angles
   - Match branch orientation to MIMICS PDF types
   - Generate parameter combinations for sensitivity analysis
5. Export combined parameters to CSV

Authors: Wanxin Yang and Russell Thomas
Last update date: 2025-11-05
================================================================================
"""

# ============================================================================
# IMPORTS
# ============================================================================
import pandas as pd
import numpy as np
from scipy import stats
import os
import argparse
import logging
from pathlib import Path
from pandarallel import pandarallel
import multiprocessing
from datetime import datetime

# ============================================================================
# MIMICS PDF FUNCTIONS
# ============================================================================

def mimics_pdf(theta_deg, pdf_type):
    """
    Calculate MIMICS probability density functions for branch orientation.
    
    Args:
        theta_deg: Angles in degrees (numpy array or float)
        pdf_type: Integer code (1-18) specifying PDF type
        
    Returns:
        PDF values at given angles
        
    PDF Types:
        1:  Uniform distribution
        2:  (sin(2θ))^4
        3-10: (sin(θ))^n where n=2,3,4,5,6,7,8,9
        11: (sin(θ-30°))^9
        12: (sin(θ+60°))^9
        13: (sin(θ))^20
        14: (sin(θ+30°))^20
        17: (sin(θ+64°))^16
        18: (sin(θ+41°))^8
    """
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


# ============================================================================
# DATA LOADING FUNCTIONS
# ============================================================================

def load_tree_info_file(info_file_path):
    """
    Load and parse RCT tree info file from RCT extraction.
    
    The info file format contains tree segments with:
    - Segment geometry (x, y, z, diameter, length)
    - Branch hierarchy (branch_order, parent_id, children)
    - Position in branch structure
    
    Args:
        info_file_path: Path to _info.txt file
        
    Returns:
        DataFrame with segment data, or None if file cannot be loaded
    """
    try:
        line_list = []
        tree_ids = []
        tree_id = 0

        with open(info_file_path, "r") as file:
            lines = file.readlines()
            for line in lines:
                # Skip comment lines (starting with #)
                if line.strip().startswith('#'):
                    continue
                
                data = line.split(", ")
                for row in data:
                    section_data = row.strip().split(", ")
                    cell_data = section_data[0].strip().split(",")
                    if len(cell_data) == 7 and all(
                        x.replace(".", "", 1).isdigit() for x in cell_data
                    ):
                        tree_id += 1
                    if len(cell_data) > 7:
                        if tree_id != 0:
                            tree_ids.append(tree_id)
                        line_list.append(cell_data)
        
        if len(line_list) < 2:
            return None
            
        df = pd.DataFrame(line_list[1:], columns=line_list[0]).astype(float)
        df.insert(0, "tree_id", tree_ids)
        # Remove row where parent_id is -1.0 (root segment)
        df = df[df["parent_id"] != -1.0]
        return df
    except Exception as e:
        print(f"Error loading {info_file_path}: {e}")
        return None


# ============================================================================
# BRANCH GEOMETRY CALCULATIONS
# ============================================================================

def calculate_branch_angles(df):
    """
    Calculate branch direction angles from segment data.
    
    For each complete branch, computes the overall orientation angle from
    first to last segment relative to vertical (z-axis).
    
    Args:
        df: DataFrame with tree segments containing branch, pos_in_branch, x, y, z
        
    Returns:
        DataFrame with added 'branch_angle_from_vertical' column (degrees)
    """
    df = df.copy()
    df['branch_angle_from_vertical'] = np.nan
   
    # Group by each complete branch
    for branch_id, group in df.groupby('branch'):
        # Sort by position in branch
        group = group.sort_values('pos_in_branch')

        if len(group) >= 2:
            # Get first and last segments
            start_segment = group.iloc[0]
            end_segment = group.iloc[-1]
           
            # Calculate overall branch direction vector
            dx = end_segment['x'] - start_segment['x']
            dy = end_segment['y'] - start_segment['y']
            dz = end_segment['z'] - start_segment['z']
           
            # Calculate length
            length = np.sqrt(dx**2 + dy**2 + dz**2)
           
            if length > 0:    
                # Branch angle from vertical
                vertical_angle = np.arccos(abs(dz) / length) * 180 / np.pi
               
                # Apply angles to all segments in this branch
                for idx in group.index:
                    df.loc[idx, 'branch_angle_from_vertical'] = vertical_angle
   
    return df


# ============================================================================
# PDF MATCHING FUNCTIONS
# ============================================================================

def calculate_pdf_match_for_tree(df_tree, census_species, min_branches=5):
    """
    Find best-matching MIMICS PDF type for each branch order in a tree.
    
    Uses volume-weighted kernel density estimation (KDE) of observed branch angles,
    then finds the MIMICS PDF with minimum KL divergence.
    
    Args:
        df_tree: DataFrame with tree segments and branch angles
        census_species: Species name (for logging)
        min_branches: Minimum number of branches required for PDF matching
        
    Returns:
        Dictionary mapping branch_order -> best_pdf_type
        e.g., {1: 4, 2: 5, 3: 6, 4: 7}
    """
    mimics_types = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 17, 18]
    theta = np.linspace(0, 180, 1000)
    
    pdf_results = {}
    
    # Group by branch_order
    for branch_order, group in df_tree.groupby('branch_order'):
        # Skip trunk (order 0) or if not enough data
        if branch_order == 0 or len(group) < 10:
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
        
        if len(angles) < min_branches:
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
            data_density = data_density / np.trapz(data_density, dx=theta[1]-theta[0])
        else:
            continue
        
        best_kl = np.inf
        best_pdf = None
        
        # Test each PDF type
        for pdf_type in mimics_types:
            mimics_pdf_values = mimics_pdf(theta, pdf_type)
            mimics_pdf_values = np.maximum(mimics_pdf_values, 0)
            integral = np.trapz(mimics_pdf_values, dx=theta[1]-theta[0])
            if integral > 0:
                mimics_pdf_normalized = mimics_pdf_values / integral
            else:
                continue
            
            # KL divergence
            data_safe = np.maximum(data_density, 1e-10)
            mimics_safe = np.maximum(mimics_pdf_normalized, 1e-10)
            kl_div = np.trapz(data_safe * np.log(data_safe / mimics_safe), dx=theta[1]-theta[0])
            
            if np.isfinite(kl_div) and kl_div < best_kl:
                best_kl = kl_div
                best_pdf = pdf_type
        
        if best_pdf is not None:
            pdf_results[branch_order] = best_pdf
    
    return pdf_results


# ============================================================================
# MAIN TREE PROCESSING FUNCTION
# ============================================================================

def process_single_tree(row, tree_attr_df, wood_density_df, 
                       rct_extraction_dir, min_branch_order=4, max_branch_order=4,
                       frequency_values=[0.43, 1.2], 
                       canopy_density_values=[0.015, 0.06, 0.24, 0.48, 0.72],
                       soil_moisture_values=[0.25, 0.5, 0.75],
                       use_matched_stems=True):
    """
    Process a single tree and generate MIMICS input parameters.
    
    This is the core processing function that:
    1. Loads RCT segment data for the tree
    2. Calculates tree-level metrics (height, volume, etc.)
    3. Computes branch statistics by order
    4. Matches branch orientations to MIMICS PDFs
    5. Generates parameter combinations for sensitivity analysis
    
    Args:
        row: Tree record from selection DataFrame
        tree_attr_df: DataFrame with tree-level attributes
        wood_density_df: DataFrame with wood density by species
        rct_extraction_dir: Path to RCT extraction directory
        min_branch_order: Minimum branch order required
        max_branch_order: Maximum branch order to process
        frequency_values: Radar frequencies for sensitivity analysis
        canopy_density_values: Canopy densities for sensitivity analysis
        soil_moisture_values: Soil moisture values for sensitivity analysis
        use_matched_stems: Whether input is from matched stems (True) or tree-attr (False)
    
    Returns:
        Tuple of (parameter_list, log_messages)
        - parameter_list: List of parameter dictionaries (one per combination)
        - log_messages: List of (level, message) tuples for logging
    """
    plot_id = row['plot_id']
    
    # -----------------------------------------------------------------------
    # SECTION 1: Identify tree and data source
    # -----------------------------------------------------------------------
    # Handle different input sources
    if use_matched_stems:
        rct_tile = row['rct_tile_best']
        rct_tree_id = row['rct_tree_id_best']
        census_species = row.get('census_species', None)
    else:
        rct_tile = row['tile']
        rct_tree_id = row['tree_id']
        census_species = None
    
    log_messages = []
    log_messages.append(('DEBUG', f"Processing tree: {plot_id} {rct_tile} {rct_tree_id} ({census_species})"))
    
    # Skip if no RCT match
    if pd.isna(rct_tile) or pd.isna(rct_tree_id):
        log_messages.append(('DEBUG', f"Skipping tree {plot_id} {rct_tree_id} - no RCT match"))
        return [], log_messages
    
    # -----------------------------------------------------------------------
    # SECTION 2: Load RCT segment data
    # -----------------------------------------------------------------------
    # Find the info file
    # Format: {country}_{plot}_raycloud_{tile_x}_{tile_y}_trees_info.txt
    # Tile format: "Tile_0_2" -> "0_2"
    tile_coords = rct_tile.replace('Tile_', '')
    # Convert plot_id to lowercase and remove hyphens for filename (e.g., "LPG-01" -> "lpg01")
    plot_lower = plot_id.lower().replace('-', '')
    
    # Get country from tree_attr_df (avoid hardcoding)
    country_matches = tree_attr_df[tree_attr_df['plot_id'] == plot_id]['country'].unique()
    if len(country_matches) > 0:
        country = country_matches[0].lower()
    else:
        country = 'unknown'  # Fallback if country not found
    
    info_file_pattern = f"{country}_{plot_lower}_raycloud_{tile_coords}_trees_info.txt"
    info_file_path = Path(rct_extraction_dir) / info_file_pattern
    
    if not info_file_path.exists():
        log_messages.append(('WARNING', f"Info file not found: {info_file_path}"))
        return [], log_messages
    
    # Load tree segments from info file
    df_all_trees = load_tree_info_file(info_file_path)
    if df_all_trees is None:
        log_messages.append(('ERROR', f"Failed to load info file: {info_file_path}"))
        return [], log_messages
    
    # Filter for this specific tree (trees are separated by blank lines in the file)
    # The info file contains all trees - we need to parse by tree structure
    # For now, let's read individual tree files
    tree_info_file = f"{country}_{plot_lower}_raycloud_{tile_coords}_trees_{int(rct_tree_id)}_info.txt"
    tree_info_path = Path(rct_extraction_dir) / f"{country}_{plot_lower}_raycloud_{tile_coords}_treesplit" / tree_info_file
    
    if not tree_info_path.exists():
        log_messages.append(('WARNING', f"Tree info file not found: {tree_info_path}"))
        return [], log_messages
    
    df_tree = load_tree_info_file(tree_info_path)
    if df_tree is None or len(df_tree) == 0:
        log_messages.append(('ERROR', f"Failed to load tree data from: {tree_info_path}"))
        return [], log_messages
    
    log_messages.append(('DEBUG', f"Loaded tree data: {len(df_tree)} segments, max branch order: {df_tree['branch_order'].max()}"))
    
    # -----------------------------------------------------------------------
    # SECTION 3: Calculate basic geometric properties
    # -----------------------------------------------------------------------
    # Calculate geometric properties for each segment
    df_tree['volume'] = np.pi * (df_tree['diameter'] / 2) ** 2 * df_tree['segment_length']
    df_tree['surface_area'] = np.pi * df_tree['diameter'] * df_tree['segment_length']
    
    # Check if tree has required branch orders
    max_order_in_tree = df_tree['branch_order'].max()
    if max_order_in_tree < min_branch_order:
        log_messages.append(('DEBUG', f"Tree {plot_id} {rct_tile} {rct_tree_id} skipped - max branch order {max_order_in_tree} < required {min_branch_order}"))
        return [], log_messages
    
    # -----------------------------------------------------------------------
    # SECTION 4: Get tree-level attributes
    # -----------------------------------------------------------------------
    tree_attr = tree_attr_df[
        (tree_attr_df['plot_id'] == plot_id) & 
        (tree_attr_df['tile'] == rct_tile) & 
        (tree_attr_df['tree_id'] == int(rct_tree_id))
    ]
    
    if len(tree_attr) == 0:
        log_messages.append(('ERROR', f"Tree attributes not found for {plot_id} {rct_tile} {rct_tree_id}"))
        return [], log_messages
    
    # Get tree attributes
    tree_attr_row = tree_attr.iloc[0]

    # Extract metadata (country and date)
    country_value = tree_attr_row.get('country', country)  # Use extracted country as fallback
    date_value = tree_attr_row.get('date', None)

    ## entire tree-level
    # tree_H = tree_attr_row['H_cloud']
    tree_H = df_tree['z'].max() - df_tree['z'].min()  # tree height calculated from point cloud
    crown_radius = tree_attr_row['crown_radius']
    dbh_m = tree_attr_row['DBH']
    total_tree_vol = df_tree['volume'].sum()
    total_tree_sa = df_tree['surface_area'].sum()

    ## stem (tree base till the first furcation)
    furcate_nodes = df_tree.loc[(df_tree.branch_order == 0.0) & (df_tree.children > 1.0), 'pos_in_branch']
    pos = furcate_nodes.min() if not furcate_nodes.empty else np.nan
    stem = df_tree[(df_tree.branch_order == 0.0) & (df_tree.pos_in_branch <= pos)]
    stem_height = stem['z'].max() - stem['z'].min()

    ## trunk (from tree base to tip of branch with order 0)
    trunk = df_tree[df_tree['branch_order'] == 0]
    trunk_vol = trunk['volume'].sum()
    trunk_sa = trunk['surface_area'].sum()
    trunk_length = trunk['segment_length'].sum()

    ## crown (parts above first furcation excluding trunk)
    # crown = df_tree.loc[~df_tree.index.isin(trunk.index)]
    crown = df_tree[df_tree['branch_order'] > 0]
    crown_thickness = crown['z'].max() - crown['z'].min() if not crown.empty else 0.0
    crown_volume = crown['volume'].sum() if not crown.empty else 0.0
    crown_surface_area = crown['surface_area'].sum() if not crown.empty else 0.0


    # -----------------------------------------------------------------------
    # SECTION 5: Calculate branch angles and match to MIMICS PDFs
    # -----------------------------------------------------------------------
    # Calculate branch angles
    df_tree = calculate_branch_angles(df_tree)
    
    # Calculate PDF matches for this tree
    pdf_matches = calculate_pdf_match_for_tree(df_tree, census_species)
    log_messages.append(('DEBUG', f"PDF matches for tree {plot_id} {rct_tile} {rct_tree_id}: {pdf_matches}"))
    
    # Check if we have PDFs for all required branch orders
    required_orders = range(1, max_branch_order + 1)
    if not all(order in pdf_matches for order in required_orders):
        log_messages.append(('WARNING', f"Missing PDF matches for tree {plot_id} {rct_tile} {rct_tree_id} - required orders: {list(required_orders)}, found: {list(pdf_matches.keys())}"))
        return [], log_messages

    # -----------------------------------------------------------------------
    # SECTION 6: Get auxiliary data (leaf stats, wood density)
    # -----------------------------------------------------------------------
    # Get leaf stats from row (already in selected_trees_df from tree_attr)
    num_leaves = row.get('num_leaves', 0) if pd.notna(row.get('num_leaves')) else 0
    leaf_area_m2 = row.get('leaf_area_m2', 0.0) if pd.notna(row.get('leaf_area_m2')) else 0.0
    log_messages.append(('DEBUG', f"Leaf stats for tree {plot_id} {rct_tile} {rct_tree_id}: {num_leaves} leaves, {leaf_area_m2} m^2"))

    # Get wood density
    wood_density_value = 0.5  # Default
    if wood_density_df is not None and census_species and pd.notna(census_species):
        density_match = wood_density_df[wood_density_df['census_species'] == census_species]
        if not density_match.empty:
            wood_density_value = density_match.iloc[0]['wood_density_mean']
            log_messages.append(('DEBUG', f"Wood density for {census_species}: {wood_density_value}"))
        else:
            log_messages.append(('DEBUG', f"No wood density found for {census_species}, using default: {wood_density_value}"))
    
    
    # -----------------------------------------------------------------------
    # SECTION 7: Calculate branch statistics by order
    # -----------------------------------------------------------------------
    # Calculate branch statistics by order
    branch_stats = {}
    for order in range(0, max_branch_order + 1):
        order_data = df_tree[df_tree['branch_order'] == order]
        
        if len(order_data) > 0:
            # Total volume and surface area of each branch order
            branch_stats[f'branch_{order}_volume'] = order_data['volume'].sum()
            branch_stats[f'branch_{order}_surface_area'] = order_data['surface_area'].sum()
            
            # Individual branch calculations
            branch_list = []
            for branch_id, branch_group in order_data.groupby('branch'):
                branch_length = branch_group['segment_length'].sum()
                branch_volume = branch_group['volume'].sum()
                
                if branch_volume > 0:
                    weighted_diameter = np.average(branch_group['diameter'], weights=branch_group['volume'])
                else:
                    weighted_diameter = branch_group['diameter'].mean()
                
                branch_list.append({
                    'length': branch_length,
                    'diameter': weighted_diameter,
                    'volume': branch_volume
                })
            
            branch_df = pd.DataFrame(branch_list)
            
            # Volume-weighted averages
            if branch_df['volume'].sum() > 0:
                branch_stats[f'branch_{order}_length'] = np.average(branch_df['length'], weights=branch_df['volume'])
                branch_stats[f'branch_{order}_diameter'] = np.average(branch_df['diameter'], weights=branch_df['volume']) * 100  # Convert to cm
            else:
                branch_stats[f'branch_{order}_length'] = branch_df['length'].mean()
                branch_stats[f'branch_{order}_diameter'] = branch_df['diameter'].mean() * 100
            
            # MIMICS crown volume: ellipsoid geometric approximation
            mimics_crown_vol = (4/3) * np.pi * (crown_radius ** 2) * ((tree_H - stem_height) / 2)
            # Number density (for non-trunk branches)
            if order != 0:
                unique_branches = order_data['branch'].nunique()
                branch_stats[f'branch_{order}_density'] = unique_branches / mimics_crown_vol if mimics_crown_vol > 0 else 0
                branch_stats[f'branch_{order}_pdf'] = pdf_matches.get(order, None)
        else:
            # Missing this order
            branch_stats[f'branch_{order}_volume'] = None
            branch_stats[f'branch_{order}_surface_area'] = None
            branch_stats[f'branch_{order}_length'] = None
            branch_stats[f'branch_{order}_diameter'] = None
            if order != 0:
                branch_stats[f'branch_{order}_density'] = None
                branch_stats[f'branch_{order}_pdf'] = None
    
    # Check if any required branch stats are None
    if any(branch_stats.get(f'branch_{order}_pdf') is None for order in required_orders):
        log_messages.append(('WARNING', f"Missing branch stats for tree {plot_id} {rct_tile} {rct_tree_id}"))
        return [], log_messages
    
    log_messages.append(('DEBUG', f"Branch stats calculated for tree {plot_id} {rct_tile} {rct_tree_id}"))
    
    # -----------------------------------------------------------------------
    # SECTION 8: Build base parameter dictionary
    # -----------------------------------------------------------------------
    # Create base parameter dictionary
    base_params = {
        'country': country_value,
        'date': date_value,
        'plot_id': plot_id,
        'tile_coords': rct_tile,
        'tree_id': int(rct_tree_id),
        'height': round(tree_H, 3),
        'total_tree_vol_m3': round(total_tree_vol, 3),
        'total_tree_sa_m2': round(total_tree_sa, 3),
        'trunk_vol_m3': round(trunk_vol, 3),
        'trunk_sa_m2': round(trunk_sa, 3),
        'trunk_diameter': round(dbh_m * 100, 3),  # DBH in cm
        'trunk_length': round(trunk_length, 3),
        'crown_vol_m3': round(crown_volume, 3),
        'crown_sa_m2': round(crown_surface_area, 3),
        'crown_radius': round(crown_radius, 3),
        'crown_thickness': round(crown_thickness, 3),
        'cstart': round(stem_height, 3),
        'leaf_density': round(leaf_area_m2 / crown_volume, 3) if crown_volume > 0 else 0,
        'census_species': census_species,
        'wood_density': round(wood_density_value, 3),
        'angle': 30,
        'rms_height': 1.5,
        'correlation_length': 17.5,
        'percent_sand': 40,
        'percent_clay': 10,
        'trunk_moisture': 0.5,
        'trunk_dry_density': round(wood_density_value, 3),
    }
    
    # -----------------------------------------------------------------------
    # SECTION 9: Add branch-specific parameters
    # -----------------------------------------------------------------------
    # Add branch parameters
    for order in range(1, max_branch_order + 1):
        base_params[f'branch_{order}_moisture'] = 0.5
        base_params[f'branch_{order}_dry_density'] = round(wood_density_value, 3)
        base_params[f'branch_{order}_length'] = round(branch_stats.get(f'branch_{order}_length', 0), 3)
        base_params[f'branch_{order}_diameter'] = round(branch_stats.get(f'branch_{order}_diameter', 0), 3)
        base_params[f'branch_{order}_density'] = round(branch_stats.get(f'branch_{order}_density', 0), 4)
        # Ensure PDF is an integer code (default to 2 if missing)
        pdf_val = branch_stats.get(f'branch_{order}_pdf', None)
        try:
            base_params[f'branch_{order}_pdf'] = int(pdf_val) if pd.notna(pdf_val) and pdf_val is not None else 2
        except Exception:
            base_params[f'branch_{order}_pdf'] = 2

    # Ensure branches 1..5 exist in the output (run_model expects up to branch_5)
    max_output_order = max(5, max_branch_order)
    for order in range(1, 6):
        if f'branch_{order}_moisture' not in base_params:
            base_params[f'branch_{order}_moisture'] = 0.5
        if f'branch_{order}_dry_density' not in base_params:
            base_params[f'branch_{order}_dry_density'] = round(wood_density_value, 3)
        if f'branch_{order}_length' not in base_params:
            base_params[f'branch_{order}_length'] = 0
        if f'branch_{order}_diameter' not in base_params:
            base_params[f'branch_{order}_diameter'] = 0
        if f'branch_{order}_density' not in base_params:
            base_params[f'branch_{order}_density'] = 0
        if f'branch_{order}_pdf' not in base_params:
            base_params[f'branch_{order}_pdf'] = 2


    # -----------------------------------------------------------------------
    # SECTION 10: Calculate volume and surface area ratios
    # -----------------------------------------------------------------------
    # Ratios (avoid division by zero)
    if total_tree_vol > 0:
        base_params['tree_sa_to_volume_ratio'] = round(total_tree_sa / total_tree_vol, 3)
        base_params['trunk_to_tree_volume_ratio'] = round(trunk_vol / total_tree_vol, 3)
        base_params['canopy_to_tree_volume_ratio'] = round(crown_volume / total_tree_vol, 3)
    else:
        base_params['tree_sa_to_volume_ratio'] = None
        base_params['trunk_to_tree_volume_ratio'] = None
        base_params['canopy_to_tree_volume_ratio'] = None

    # Per-order ratios for orders 0 till max_output_order
    for order in range(0, max_output_order + 1):
        order_vol = branch_stats.get(f'branch_{order}_volume', 0) or 0
        order_sa = branch_stats.get(f'branch_{order}_surface_area', 0) or 0
        if total_tree_vol > 0:
            base_params[f'order_{order}_to_tree_volume_ratio'] = round(order_vol / total_tree_vol, 3)
        else:
            base_params[f'order_{order}_to_tree_volume_ratio'] = None

        if order_vol > 0:
            base_params[f'order_{order}_sa_to_volume_ratio'] = round(order_sa / order_vol, 3) if order_sa is not None else None
        else:
            base_params[f'order_{order}_sa_to_volume_ratio'] = None

    # Ensure the PDF fields are integers (final safety)
    for order in range(1, max_output_order + 1):
        key = f'branch_{order}_pdf'
        try:
            base_params[key] = int(base_params.get(key, 2))
        except Exception:
            base_params[key] = 2

    # -----------------------------------------------------------------------
    # SECTION 11: Generate parameter combinations for sensitivity analysis
    # -----------------------------------------------------------------------
    # Generate parameter combinations for sensitivity analysis
    # Each tree gets combinations based on input parameter lists
    results = []
    for freq_val in frequency_values:
        for canopy_dens_val in canopy_density_values:
            for soil_moist_val in soil_moisture_values:
                params = base_params.copy()
                params['frequency'] = freq_val
                params['canopy_density'] = round(canopy_dens_val, 5)
                params['soil_moisture'] = soil_moist_val
                results.append(params)

    log_messages.append(('DEBUG', f"Generated {len(results)} parameter combinations for tree {plot_id} {rct_tile} {rct_tree_id}"))

    return results, log_messages


# ============================================================================
# COMMAND LINE INTERFACE AND MAIN WORKFLOW
# ============================================================================

def main():
    """
    Main entry point for MIMICS parameter generation.
    
    Workflow:
    1. Parse command line arguments
    2. Load and validate input data (single plot only)
    3. Filter trees by DBH and height thresholds
    4. Process trees in parallel to generate parameters
    5. Save combined results to CSV
    """
    
    # ========================================================================
    # STEP 1: Parse command line arguments
    # ========================================================================
    parser = argparse.ArgumentParser(
        description='Generate MIMICS model input parameters for RayCloudTools extracted trees',
        epilog='''
Example usage:
  # Process trees using census matched stems:
  python %(prog)s -t tree_attr.csv -r ../rct_extraction/ -c census_matched_stems.csv -w wood_density.csv -odir ../mimics_inputs/ 
  
  # Process all in-plot trees without census matching:
  python %(prog)s -t tree_attr.csv -r ../rct_extraction/ -odir ../mimics_inputs/ 

  # With custom thresholds and sensitivity analysis parameters:
  python %(prog)s -t tree_attr.csv -r ../rct_extraction/ -d 10.0 -H 3.0 -f 0.43 1.2 5.4 -s 0.25 0.5 0.75 --canopy_density 0.015 0.06 0.24
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument(
        '-t', '--tree_attr_csv',
        type=str,
        required=True,
        help='Path to tree-level attributes CSV file'
    )
    
    parser.add_argument(
        '-r', '--rct_extraction_dir',
        type=str,
        required=True,
        help='Directory containing RCT extraction results'
    )
    
    parser.add_argument(
        '-odir', '--output_dir',
        type=str,
        default='.',
        help='Directory to store the output MIMICS input parameters CSV file (default: current directory)'
    )
    
    # Optional arguments
    parser.add_argument(
        '-c', '--census_matched_stems_csv',
        type=str,
        default=None,
        help='Path to matched stems CSV file for this plot (optional). If not provided, all in-plot trees from tree_attr_csv will be processed.'
    )
    
    parser.add_argument(
        '-w', '--wood_density_csv',
        type=str,
        default=None,
        help='Path to wood density CSV file (optional)'
    )
    
    parser.add_argument(
        '-a', '--area',
        type=float,
        default=1.0,
        help='Plot area in hectares (used for canopy density calculation if not provided)'
    )
    
    parser.add_argument(
        '--min_branch_order',
        type=int,
        default=3,
        help='Minimum branch order required for tree inclusion (default: 3)'
    )
    
    parser.add_argument(
        '--max_branch_order',
        type=int,
        default=4,
        help='Maximum branch order to process in MIMICS (default: 4)'
    )
    
    parser.add_argument(
        '-d', '--DBH_threshold',
        type=float,
        default=10.0,
        help='Minimum DBH threshold in cm for tree inclusion (default: 10 cm)'
    )
    
    parser.add_argument(
        '-H', '--H_threshold',
        type=float,
        default=3.0,
        help='Minimum height threshold in meters for tree inclusion (default: 3 m)'
    )
    
    parser.add_argument(
        '-n', '--n_cores',
        type=int,
        default=max(1, multiprocessing.cpu_count() - 1),
        help='Number of CPU cores for parallel processing (default: max available - 1)'
    )
    
    parser.add_argument(
        '-f', '--radar_frequency',
        type=float,
        nargs='*',
        default=[0.43, 1.2, 5.4],
        help='Radar operating frequency for sensitivity analysis (space-separated list or single value), default: 0.43 1.2 5.4'
    )
    
    parser.add_argument(
        '--canopy_density',
        type=float,
        nargs='*',
        default=None,
        help='Canopy density (number of trees per square meter) (space-separated list or single value). If not provided, will be calculated from selected trees count / plot area.'
    )
    
    parser.add_argument(
        '-s', '--soil_moisture',
        type=float,
        nargs='*',
        default=[0.25, 0.5, 0.75],
        help='Soil moisture values for sensitivity analysis (space-separated list or single value), default: 0.25 0.5 0.75'
    )
    
    args = parser.parse_args()
    
    # ========================================================================
    # STEP 2: Load input data and validate plot consistency
    # ========================================================================
    # Set up logging
    output_dir = args.output_dir if args.output_dir else '.'
    
    # Load tree attributes first to determine and validate plot_id
    tree_attr_df = pd.read_csv(args.tree_attr_csv)
    
    # Get unique plot_ids from tree_attr_csv
    tree_attr_plot_ids = tree_attr_df['plot_id'].unique()
    
    if len(tree_attr_plot_ids) == 0:
        print(f"ERROR: No plot_id found in tree_attr_csv file!")
        print(f"  File: {args.tree_attr_csv}")
        raise ValueError("tree_attr_csv file contains no plot_id data")
    
    if len(tree_attr_plot_ids) > 1:
        print(f"ERROR: tree_attr_csv file contains multiple plot_ids: {tree_attr_plot_ids}")
        print(f"  File: {args.tree_attr_csv}")
        print(f"  This script processes only one plot at a time.")
        print(f"  Please filter the tree_attr_csv file to contain only one plot.")
        raise ValueError(f"tree_attr_csv contains multiple plots: {tree_attr_plot_ids}")
    
    # Single plot_id from tree_attr
    plot_id = tree_attr_plot_ids[0]
    
    # If matched stems file is provided, validate it has the same plot_id
    if args.census_matched_stems_csv:
        matched_stems_df = pd.read_csv(args.census_matched_stems_csv)
        matched_plot_ids = matched_stems_df['plot_id'].unique()
        
        if len(matched_plot_ids) == 0:
            print(f"ERROR: No plot_id found in census_matched_stems_csv file!")
            print(f"  File: {args.census_matched_stems_csv}")
            raise ValueError("census_matched_stems_csv file contains no plot_id data")
        
        if len(matched_plot_ids) > 1:
            print(f"ERROR: census_matched_stems_csv file contains multiple plot_ids: {matched_plot_ids}")
            print(f"  File: {args.census_matched_stems_csv}")
            print(f"  This script processes only one plot at a time.")
            print(f"  Please filter the census_matched_stems_csv file to contain only one plot.")
            raise ValueError(f"census_matched_stems_csv contains multiple plots: {matched_plot_ids}")
        
        matched_plot_id = matched_plot_ids[0]
        
        # Check if plot_ids match
        if matched_plot_id != plot_id:
            print(f"ERROR: Plot ID mismatch between input files!")
            print(f"  tree_attr_csv plot_id: {plot_id}")
            print(f"  census_matched_stems_csv plot_id: {matched_plot_id}")
            print(f"  Both files must be from the same plot.")
            print(f"  Please check your input files and rerun.")
            raise ValueError(f"Plot ID mismatch: tree_attr has '{plot_id}' but matched_stems has '{matched_plot_id}'")
    
    # Extract country from tree_attr_df (assuming it has a 'country' column)
    country = tree_attr_df['country'].iloc[0] if 'country' in tree_attr_df.columns else 'unknown'
    
    # Extract date from tree_attr_df if available
    has_date = 'date' in tree_attr_df.columns
    date_value = None
    if has_date:
        date_value = tree_attr_df['date'].iloc[0]
        if pd.isna(date_value):
            has_date = False
    
    # ========================================================================
    # STEP 3: Set up logging with timestamp
    # ========================================================================
    # Determine output filename based on date availability (needed for log file naming)
    if has_date:
        output_filename = f'mimics_inputs_{date_value}_{country}_{plot_id}.csv'
    else:
        output_filename = f'mimics_inputs_{country}_{plot_id}.csv'
    
    # Create log filename from output filename with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = output_filename.replace('.csv', f'_{timestamp}.log')
    log_file = os.path.join(output_dir, log_filename)
    
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()  # Also log to console
        ]
    )
    
    logger = logging.getLogger(__name__)
    
    print("\n" + "=" * 80)
    print("MIMICS PARAMETER GENERATION")
    print("=" * 80)
    print(f"Plot: {plot_id} ({country})")
    print(f"Log file: {log_filename}")
    print("=" * 80)
    
    logger.info("=" * 80)
    logger.info("MIMICS PARAMETER GENERATION")
    logger.info("=" * 80)
    logger.info(f"Plot: {plot_id} ({country})")
    
    # ========================================================================
    # STEP 4: Log configuration and input arguments
    # ========================================================================
    # Log configuration
    mode = "census-matched stems" if args.census_matched_stems_csv else "all in-plot trees"
    
    print(f"\n[1/5] Configuration")
    print(f"  Mode: {mode}")
    print(f"  Tree attributes: {args.tree_attr_csv}")
    if args.census_matched_stems_csv:
        print(f"  Census matches: {args.census_matched_stems_csv}")
    print(f"  RCT extraction: {args.rct_extraction_dir}")
    if args.wood_density_csv:
        print(f"  Wood density: {args.wood_density_csv}")
    print(f"  Filters: DBH ≥ {args.DBH_threshold} cm, Height ≥ {args.H_threshold} m")
    print(f"  Branch orders: {args.min_branch_order}-{args.max_branch_order}")
    print(f"  Output: {args.output_dir}")
    print(f"  Cores: {args.n_cores}")
    
    logger.info(f"Configuration:")
    logger.info(f"  Mode: {mode}")
    logger.info(f"  Tree attributes: {args.tree_attr_csv}")
    if args.census_matched_stems_csv:
        logger.info(f"  Census matches: {args.census_matched_stems_csv}")
    logger.info(f"  RCT extraction: {args.rct_extraction_dir}")
    if args.wood_density_csv:
        logger.info(f"  Wood density: {args.wood_density_csv}")
    logger.info(f"  Filters: DBH ≥ {args.DBH_threshold} cm, Height ≥ {args.H_threshold} m")
    logger.info(f"  Branch orders: {args.min_branch_order}-{args.max_branch_order}")
    logger.info(f"  Radar frequencies: {args.radar_frequency} GHz")
    logger.info(f"  Soil moisture values: {args.soil_moisture}")
    if args.canopy_density:
        logger.info(f"  Canopy density values: {args.canopy_density} trees/m²")
    logger.info(f"  Output: {args.output_dir}")
    logger.info(f"  Cores: {args.n_cores}")
    
    # ========================================================================
    # STEP 5: Load additional data files
    # ========================================================================
    print(f"\n[2/5] Loading data files...")
    
    logger.info("Loading data files...")
    logger.info(f"  Tree attributes: {len(tree_attr_df)} rows")
    
    # Optionally load matched stems
    use_matched_stems = args.census_matched_stems_csv is not None
    if use_matched_stems:
        # matched_stems_df already loaded during validation
        print(f"  ✓ Tree attributes: {len(tree_attr_df)} rows")
        print(f"  ✓ Census matches: {len(matched_stems_df)} rows")
        logger.info(f"  Census matches: {len(matched_stems_df)} rows")
    else:
        matched_stems_df = None
        print(f"  ✓ Tree attributes: {len(tree_attr_df)} rows")
        logger.info(f"  Mode: Using all in-plot trees")

    wood_density_df = pd.read_csv(args.wood_density_csv) if args.wood_density_csv else None
    if wood_density_df is not None:
        print(f"  ✓ Wood density: {len(wood_density_df)} species")
        logger.info(f"  Wood density: {len(wood_density_df)} species")
    else:
        print(f"  ✓ Wood density: Using default (0.5 g/cm³)")
        logger.info(f"  Wood density: Using default (0.5 g/cm³)")
    
    # ========================================================================
    # STEP 6: Filter and select trees to process
    # ========================================================================
    print(f"\n[3/5] Selecting trees...")
    
    logger.info(f"Selecting trees for plot {plot_id}...")
    
    if use_matched_stems:
        # Use matched stems approach
        logger.info(f"  Using census matched stems")
        
        # Filter for this plot and only those with RCT matches (plot_id already validated)
        matched_subset = matched_stems_df[
            (matched_stems_df['rct_tile_best'].notna()) &
            (matched_stems_df['rct_tree_id_best'].notna())
        ].copy()
        
        logger.info(f"  Matched stems with RCT: {len(matched_subset)}")
        
        # Deduplicate: keep only the first (best) match per RCT tree
        # This prevents one RCT tree from being processed multiple times with different census species
        matched_subset_dedup = matched_subset.drop_duplicates(
            subset=['plot_id', 'rct_tile_best', 'rct_tree_id_best'], 
            keep='first'
        )
        
        duplicates_removed = len(matched_subset) - len(matched_subset_dedup)
        if duplicates_removed > 0:
            logger.info(f"  Removed {duplicates_removed} duplicate matches")
        
        # Extract (plot_id, tile, tree_id) pairs from deduplicated matched stems
        tree_keys = matched_subset_dedup[['plot_id', 'rct_tile_best', 'rct_tree_id_best']].rename(
            columns={'rct_tile_best': 'tile', 'rct_tree_id_best': 'tree_id'}
        )
        
        # Find corresponding trees in tree_attr_df
        selected_trees_df = tree_attr_df.merge(
            tree_keys,
            on=['plot_id', 'tile', 'tree_id'],
            how='inner'
        )
        
        # Add census_species and rct_tile_best, rct_tree_id_best from deduplicated matched_stems
        selected_trees_df = selected_trees_df.merge(
            matched_subset_dedup[['plot_id', 'rct_tile_best', 'rct_tree_id_best', 'census_species']],
            left_on=['plot_id', 'tile', 'tree_id'],
            right_on=['plot_id', 'rct_tile_best', 'rct_tree_id_best'],
            how='left'
        )
        # Keep rct_tile_best and rct_tree_id_best for process_single_tree function
        
        logger.info(f"  Found in tree_attr: {len(selected_trees_df)}")
    else:
        # Use tree-attr approach
        logger.info(f"  Using all in-plot trees from tree-attr")
        
        # Filter only in-plot trees (plot_id already validated to be single value)
        selected_trees_df = tree_attr_df[tree_attr_df['in_plot'] == True].copy()
        
        logger.info(f"  In-plot trees: {len(selected_trees_df)}")
    
    # Apply DBH and height thresholds
    # Convert DBH from m to cm if needed
    if 'DBH' in selected_trees_df.columns:
        # Check if DBH is in meters (values < 1) or cm (values > 10)
        if selected_trees_df['DBH'].mean() < 1:
            selected_trees_df['DBH_cm'] = selected_trees_df['DBH'] * 100
        else:
            selected_trees_df['DBH_cm'] = selected_trees_df['DBH']
    else:
        print("  ⚠ WARNING: DBH column not found in data")
        logger.warning("DBH column not found in data")
        selected_trees_df['DBH_cm'] = 0
    
    # Filter by thresholds
    before_filter = len(selected_trees_df)
    selected_trees_df = selected_trees_df[
        (selected_trees_df['DBH_cm'] >= args.DBH_threshold) &
        (selected_trees_df['height'] >= args.H_threshold)
    ]
    after_filter = len(selected_trees_df)
    
    print(f"  Candidates: {before_filter} trees")
    print(f"  After filtering (DBH ≥ {args.DBH_threshold} cm, H ≥ {args.H_threshold} m): {after_filter} trees")
    
    logger.info(f"  Before filtering: {before_filter} trees")
    logger.info(f"  After DBH/Height filtering: {after_filter} trees ({before_filter - after_filter} filtered out)")
    
    if len(selected_trees_df) == 0:
        error_msg = f"No trees remaining after filtering!"
        print(f"\n  ✗ ERROR: {error_msg}")
        print(f"    Check: DBH ≥ {args.DBH_threshold} cm and Height ≥ {args.H_threshold} m")
        logger.error(error_msg)
        logger.error(f"  Thresholds: DBH ≥ {args.DBH_threshold} cm, Height ≥ {args.H_threshold} m")
        raise ValueError(f"No trees remaining for plot {plot_id} after applying DBH and height thresholds")
    
    # ========================================================================
    # STEP 7: Calculate canopy density and parameter combinations
    # ========================================================================
    # Calculate canopy density for this plot if not provided
    if args.canopy_density is None:
        # Calculate from tree count / area (convert ha to m^2)
        plot_area_m2 = args.area * 10000
        calculated_canopy_density = len(selected_trees_df) / plot_area_m2
        canopy_density_values = [calculated_canopy_density]
        print(f"  Canopy density: {calculated_canopy_density:.6f} trees/m² (auto-calculated)")
        logger.info(f"  Canopy density: {calculated_canopy_density:.6f} trees/m² ({len(selected_trees_df)} trees / {args.area} ha)")
    else:
        canopy_density_values = args.canopy_density
        print(f"  Canopy density: {canopy_density_values} trees/m² (user-specified)")
        logger.info(f"  Canopy density values: {canopy_density_values} trees/m²")
    
    # Calculate expected combinations
    expected_combinations = len(args.radar_frequency) * len(canopy_density_values) * len(args.soil_moisture)
    print(f"  Parameter sets per tree: {expected_combinations} ({len(args.radar_frequency)} freq × {len(canopy_density_values)} density × {len(args.soil_moisture)} moisture)")
    logger.info(f"  Expected combinations per tree: {expected_combinations}")
    
    # ========================================================================
    # STEP 8: Process trees in parallel
    # ========================================================================
    print(f"\n[4/5] Processing {len(selected_trees_df)} trees (using {args.n_cores} cores)...")
    
    logger.info(f"Processing trees...")
    logger.info(f"  Trees to process: {len(selected_trees_df)}")
    logger.info(f"  CPU cores: {args.n_cores}")
    
    pandarallel.initialize(nb_workers=args.n_cores, progress_bar=True, verbose=0)
    
    # Process trees in parallel
    results = selected_trees_df.parallel_apply(
        lambda row: process_single_tree(
            row, 
            tree_attr_df, 
            wood_density_df,
            args.rct_extraction_dir,
            args.min_branch_order,
            args.max_branch_order,
            args.radar_frequency,
            canopy_density_values,
            args.soil_moisture,
            use_matched_stems
        ),
        axis=1
    )
    
    # ========================================================================
    # STEP 9: Collect and log results
    # ========================================================================
    print(f"\n  Collecting results...")
    
    # Flatten results (each tree returns a tuple of parameter sets and log messages)
    all_params = []
    successful_trees = 0
    failed_trees = 0
    all_log_messages = []
    
    for result_tuple in results:
        result_list, log_messages = result_tuple
        all_log_messages.extend(log_messages)
        
        if result_list and len(result_list) > 0:
            all_params.extend(result_list)
            successful_trees += 1
        else:
            failed_trees += 1
    
    # Log all collected messages (only warnings and errors for cleaner output)
    for level, message in all_log_messages:
        if level in ['WARNING', 'ERROR']:
            if level == 'WARNING':
                logger.warning(message)
            elif level == 'ERROR':
                logger.error(message)
    
    print(f"\n  Results:")
    print(f"    ✓ Successful: {successful_trees} trees")
    if failed_trees > 0:
        print(f"    ✗ Failed/skipped: {failed_trees} trees")
    print(f"    → Parameter sets: {len(all_params)} (expected: {successful_trees * expected_combinations})")
    
    logger.info(f"Processing complete:")
    logger.info(f"  Successful: {successful_trees} trees")
    logger.info(f"  Failed/skipped: {failed_trees} trees")
    logger.info(f"  Parameter sets: {len(all_params)} (expected: {successful_trees * expected_combinations})")
    
    if len(all_params) == 0:
        error_msg = "No valid parameter sets generated!"
        print(f"\n  ✗ ERROR: {error_msg}")
        print(f"    Possible reasons:")
        print(f"      - Trees don't meet min branch order requirement ({args.min_branch_order})")
        print(f"      - RCT info files not found")
        print(f"      - Insufficient branches for PDF matching")
        logger.error(error_msg)
        logger.error(f"  Check: min_branch_order={args.min_branch_order}, RCT extraction directory")
        raise ValueError(f"No valid MIMICS parameters generated for plot {plot_id}")
    
    # ========================================================================
    # STEP 10: Save output CSV file
    # ========================================================================
    print(f"\n[5/5] Saving results...")
    
    logger.info("Saving output...")
    
    df_output = pd.DataFrame(all_params)
    
    print(f"  Output: {df_output.shape[0]} rows × {df_output.shape[1]} columns")
    logger.info(f"  Output shape: {df_output.shape[0]} rows × {df_output.shape[1]} columns")
    
    # Output filename already determined earlier (for log file naming)
    output_file = os.path.join(args.output_dir, output_filename)
    
    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        logger.info(f"  Created output directory: {args.output_dir}")
    
    df_output.to_csv(output_file, index=False)
    
    file_size_kb = os.path.getsize(output_file) / 1024
    print(f"  ✓ Saved: {output_file}")
    print(f"  Size: {file_size_kb:.1f} KB")
    
    logger.info(f"  Output file: {output_file}")
    logger.info(f"  File size: {file_size_kb:.1f} KB")
    
    print(f"\n{'='*80}")
    print(f"✓ COMPLETED SUCCESSFULLY")
    print(f"{'='*80}")
    print(f"  Plot: {plot_id}")
    print(f"  Trees processed: {successful_trees}/{len(selected_trees_df)}")
    print(f"  Parameters generated: {len(all_params)}")
    print(f"  Output: {output_file}")
    print("=" * 80 + "\n")
    
    logger.info("=" * 80)
    logger.info("COMPLETED SUCCESSFULLY")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
