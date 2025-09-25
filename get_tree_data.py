import pandas as pd
import numpy as np
from pathlib import Path
from plyfile import PlyData

def leaf_area(leaf_file):
    """
    Calculate the surface area of a leaf from a PLY mesh file.
    
    Reads a PLY file containing leaf mesh data, computes the area of each triangular face,
    and returns both the total surface area and area.
    Args:
        leaf_file: Path to PLY file containing leaf mesh data
    Returns:
        tuple: (total_area, total_area * 1000) - surface area in two units
    """
    mesh = PlyData.read(leaf_file)
    vertices = np.vstack([mesh['vertex'][axis] for axis in ['x', 'y', 'z']]).T
    faces = np.vstack(mesh['face']['vertex_indices']) if 'face' in mesh else None
    if faces is None:
        return 0, 0
    v0, v1, v2 = vertices[faces[:,0]], vertices[faces[:,1]], vertices[faces[:,2]]
    tri_areas = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
    total_area = tri_areas.sum()
    return total_area, total_area * 1000

def treeinfo_attributes_tree(tree_file):
    """
    Extracts per-tree tree attributes from a tree file generated using treeinfo and returns a DataFrame.
    Can be run on a 'forest' or single treefile.
    Args:
    tree_file (str): The path to the treefile.
    Returns:
    pandas.DataFrame: A DataFrame containing tree attributes. If the input file is a forest, the DataFrame will contain a column 'tree_id' with the tree ID.
    """
    line_list = []
    with open(tree_file, "r") as file:
        lines = file.readlines()
        line_count = 0
        for line in lines:
            data = line.split(", ")
            for row in data:
                section_data = row.strip().split(", ")
                cell_data = section_data[0].strip().split(",")
                if len(cell_data) == 7:
                    line_list.append(cell_data)
            line_count += 1
    
    if not line_list:
        return pd.DataFrame()
    
    df = pd.DataFrame(line_list[1:], columns=line_list[0]).astype(float)
    if len(df) > 1:
        df.insert(0, "tree_id", range(1, len(df) + 1))
    return df

def treeinfo_attributes_segment(tree_file):
    """
    Extracts per-segment attributes of a tree file generated using treeinfo and returns a DataFrame.
    Can be run on a 'forest' or single treefile.
    Parameters:
    tree_file (str): The path to the treefile created using treeinfo.
    Returns:
    pandas.DataFrame: A DataFrame containing segment attributes. If the input file is a forest, the DataFrame will contain a column 'tree_id' with the tree ID.
    """
    line_list = []
    tree_ids = []
    tree_id = 0
    with open(tree_file, "r") as file:
        lines = file.readlines()
        line_count = 0
        for line in lines:
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
            line_count += 1
    
    if not line_list:
        return pd.DataFrame()
    
    df = pd.DataFrame(line_list[1:], columns=line_list[0]).astype(float)
    if tree_ids:  # Only add tree_id column if we have tree_ids
        df.insert(0, "tree_id", tree_ids)
    # remove row where parent_id is -1.0
    df = df[df["parent_id"] != -1.0]
    return df

def extract_good_trees(tree_files_dir, matched_stems_file, output_dir="output"):
    '''
    Extract and process data for trees marked as 'good' in the matched stems file.
    
    Searches through tree info files to find matches for good trees, extracts segment and tree
    attributes, adds leaf area data if available, and saves combined results grouped by 
    plot_id and tile coordinates as separate CSV files.
    
    Args:
        tree_files_dir: Directory containing tree info files to search
        matched_stems_file: CSV file with tree suitability ratings
        output_dir: Output directory for combined CSV files (default: "output")
    '''
    tree_files_path = Path(tree_files_dir)
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Load the matched stems file and filter for good trees
    matched_df = pd.read_csv(matched_stems_file)
    good_trees = matched_df[(matched_df['suitability'].str.lower() == 'good')].copy()
    
    print(f"Found {len(good_trees)} good trees to process")
    
    # Find all tree info files
    # '*/' SEARCH All subdirectoires. Then '/*_trees_*_info.txt' Searches each file
    tree_files = list(tree_files_path.glob("*/*_trees_*_info.txt"))
    print(f"Found {len(tree_files)} tree files to search")
    
    all_results = []
    
    for _, good_tree in good_trees.iterrows():
        target_plot = good_tree['plot_id']
        target_tile = good_tree['rct_tile_best']
        target_tree_id = good_tree['rct_tree_id_best']
        
        print(f"Looking for {target_plot} {target_tile} tree {target_tree_id}")
        
        # Find the matching tree file
        found_file = None
        for tree_file in tree_files:
            if (f"_{target_plot}_" in tree_file.name and 
                f"_{target_tile.replace('_', '_')}_" in tree_file.name and
                f"_trees_{target_tree_id}_" in tree_file.name):
                found_file = tree_file
                break
        
        if not found_file:
            print(f"  File not found")
            continue
            
        print(f"  Found: {found_file.name}")
        
        # Get segment data using the working function
        segment_df = treeinfo_attributes_segment(str(found_file))
        if segment_df.empty:
            print(f"  No segment data")
            continue
        
        # Add data from matched stems
        for col in good_tree.index:
            segment_df[col] = good_tree[col]
        
        segment_df['file_name'] = found_file.name
        
        # Clean up column names
        segment_df = segment_df.drop('tree_id', axis=1, errors='ignore')
        segment_df = segment_df.rename(columns={
            'rct_tree_id_best': 'tree_id',
            'rct_tile_best': 'tile_coords',
            'rct_dbh_best': 'dbh_best'
        })
        
        # Get tree data using the treeinfo function
        tree_df = treeinfo_attributes_tree(str(found_file))
        if not tree_df.empty:
            print(f"  Found tree data with {len(tree_df)} rows")
            # For individual tree files, add the tree data as new columns
            for col in tree_df.columns:
                if col not in segment_df.columns:
                    segment_df[col] = tree_df[col].iloc[0] if len(tree_df) > 0 else None
        else:
            print(f"No tree data")
        
        # Get leaf data
        treesplit_dir = found_file.parent
        leaf_pattern = found_file.name.replace('_trees_', '_segmented_').replace('_info.txt', '_leaves.ply')
        leaf_files = list(treesplit_dir.glob(leaf_pattern))
        
        if leaf_files:
            try:
                total_area, num_leaves = leaf_area(leaf_files[0])
                segment_df['leaf_area_m2'] = total_area
                segment_df['num_leaves'] = num_leaves
                print(f"  Added leaf data: {total_area:.4f} m²")
            except Exception as e:
                print(f"  Error processing leaf: {e}")
        
        all_results.append(segment_df)
        print(f"  Added {len(segment_df)} segments")
    
    # Group results by plot_id and tile_coords for separate files
    plot_tile_groups = {}
    for result_df in all_results:
        if not result_df.empty:
            plot_id = result_df['plot_id'].iloc[0]
            tile_coords = result_df['tile_coords'].iloc[0]  # Updated column name
            group_key = f"angola_{plot_id}_{tile_coords}"
            
            if group_key not in plot_tile_groups:
                plot_tile_groups[group_key] = []
            plot_tile_groups[group_key].append(result_df)
    
    # Save each group as separate CSV file
    total_segments = 0
    for group_name, group_dfs in plot_tile_groups.items():
        combined_group = pd.concat(group_dfs, ignore_index=True)
        csv_file = output_path / f"{group_name}_combined.csv"
        combined_group.to_csv(csv_file, index=False)
        
        total_segments += len(combined_group)
        print(f"{group_name}: {len(combined_group)} segments saved")
    
    print(f"\nProcessed {len(all_results)} trees")
    print(f"Total segments: {total_segments}")
    print(f"Created {len(plot_tile_groups)} output files")

if __name__ == "__main__":
    tree_files_dir = "/home/ucfargt@ad.ucl.ac.uk/Documents/Thesis/mimics/TLS-QSM_results_angola_bicuar_tree/angola_results/rct_extraction"
    matched_stems_file = "/home/ucfargt@ad.ucl.ac.uk/Documents/Thesis/mimics/TLS-QSM_results_angola_bicuar_tree/dataset/matched_stems_picked.csv"
    
    extract_good_trees(tree_files_dir, matched_stems_file, "output")
