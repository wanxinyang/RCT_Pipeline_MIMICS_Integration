#%%
"""
================================================================================
RayCloudTools Tree-Level Attributes Extraction
================================================================================

This script extracts tree-level attributes from RayCloudTools (RCT) outputs,
including tree measurements, leaf area calculations, and spatial filtering.

Main functionalities:
1. Read tree attributes from RCT treeinfo and tree files
2. Merge tree data and filter duplicates across tiles
3. Calculate plot boundaries from scan positions
4. Filter in-plot trees based on spatial boundaries
5. Calculate leaf area from mesh files (parallelized)
6. Export results to CSV with comprehensive metadata

Author: Wanxin Yang
Last Modified: 2025-10-29
"""

#%%
# ============================================================================
# IMPORT LIBRARIES
# ============================================================================
# Standard library
import sys
import os
from glob import glob
from pathlib import Path
import argparse
from datetime import datetime
import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

# Third-party libraries
import numpy as np
import pandas as pd
import geopandas as gp
from shapely.geometry import Point
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from plyfile import PlyData
from tqdm import tqdm

#%%
# ============================================================================
# COMMAND-LINE ARGUMENT PARSING
# ============================================================================
def parse_args():
    parser = argparse.ArgumentParser(
        description='Extract tree-level attributes from RayCloudTools outputs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Usage example:\n"
            "python extract_rct_tree_attrs.py\n"
            "  -m /path/to/matrix_dir\n"
            "  -t /path/to/tile_index.dat\n"
            "  -r /path/to/rct_extraction_dir\n"
            "  -c Angola -p P02 -d 2024-02-23\n"
            "  -o /path/to/output_dir -f /path/to/fig_dir\n"
        )
    )
    
    # Required input paths
    parser.add_argument('-m', '--dir_matrix', type=str, required=True,
                        help='Path to directory containing SOP matrix files (*.DAT)')
    parser.add_argument('-t', '--fp_tile', type=str, required=True,
                        help='Path to tile index file (tile_index.dat)')
    parser.add_argument('-r', '--dir_rct_extraction', type=str, required=True,
                        help='Path to Raycloudtools extraction results directory')

    # Required plot metadata
    parser.add_argument('-c', '--country', type=str, required=True,
                        help='Country name (REQUIRED), e.g., Angola')
    parser.add_argument('-p', '--plotid', type=str, required=True,
                        help='Plot ID (REQUIRED), e.g., P02')
    parser.add_argument('-d', '--date', type=str, default=None,
                        help='Date in YYYY-MM-DD format (e.g. 2024-02-23, default: None)')
    
    # Optional output directories (default to current directory)
    parser.add_argument('-o', '--odir', type=str, default='.',
                        help='Output directory for tree attributes CSV files (default: current directory)')
    parser.add_argument('-f', '--figdir', type=str, default='.',
                        help='Output directory for figures (default: current directory)')
    
    # Tiling parameters
    parser.add_argument('--rct_tile_len', type=int, default=50,
                        help='Raycloudtools tile length (metres) used to generated the extraction results')
    parser.add_argument('--rct_tile_overlap', type=int, default=10,
                        help='Raycloudtools tile overlap (metres) used to generate the extraction results')

    # Processing parameters
    max_cores = max(1, multiprocessing.cpu_count() - 1)
    parser.add_argument('--n_cores', type=int, default=max_cores-1,
                        help=f'Number of parallel cores for leaf processing (default: {max_cores-1})')
    
    return parser.parse_args()

#%%
# ============================================================================
# LOGGING SETUP
# ============================================================================
def setup_logging(odir, country, plotid, date):
    """
    Set up logging to both file and console.
    
    Args:
        odir (Path): Output directory for log files
        country (str): Country name for log filename
        plotid (str): Plot ID for log filename
        date (str or None): Date string for log filename
        
    Returns:
        logging.Logger: Configured logger instance
    """
    # Create log filename with timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    date_str = date if date else 'NoDate'
    log_filename = odir / f'{date_str}_{country}_{plotid}_tree_extraction_{timestamp}.log'
    
    # Configure logging
    logger = logging.getLogger('TreeExtraction')
    logger.setLevel(logging.INFO)
    
    # Remove existing handlers to avoid duplicates
    logger.handlers.clear()
    
    # File handler
    fh = logging.FileHandler(log_filename)
    fh.setLevel(logging.INFO)
    
    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    
    # Formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s',
                                  datefmt='%Y-%m-%d %H:%M:%S')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    # Add handlers
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    logger.info("="*80)
    logger.info("RayCloudTools Tree Extraction Pipeline Started")
    logger.info("="*80)
    logger.info(f"Log file: {log_filename}")
    
    return logger

#%%
# ============================================================================
# CONFIGURATION AND DATA LOADING
# ============================================================================
# Check if running interactively or from command line
if __name__ == '__main__':
    args = parse_args()
    
    # Validate required metadata
    if not args.country:
        raise ValueError("Country information is required. Please provide --country argument.")
    if not args.plotid:
        raise ValueError("Plot ID information is required. Please provide --plotid argument.")
    
    # Extract arguments
    country = args.country
    plotid = args.plotid
    date = args.date  # Can be None
    rct_tile_len = args.rct_tile_len
    rct_tile_overlap = args.rct_tile_overlap
    dir_matrix = Path(args.dir_matrix)
    dat_files = list(dir_matrix.glob('*.DAT')) + list(dir_matrix.glob('*.dat'))
    fp_tile = Path(args.fp_tile)
    parent_dir = Path(args.dir_rct_extraction)
    odir = Path(args.odir)
    figdir = Path(args.figdir)
    num_workers = args.n_cores
    
    # Validate output directories exist
    if not odir.exists():
        raise ValueError(f"Output directory does not exist: {odir}\nPlease create it first or specify an existing directory.")
    if not figdir.exists():
        raise ValueError(f"Figure directory does not exist: {figdir}\nPlease create it first or specify an existing directory.")
    
    # Setup logging
    logger = setup_logging(odir, country, plotid, date)
    
else:
    # Default values for interactive use
    country = 'Angola'
    plotid = 'P02'
    date = '2024-02-23'
    rct_tile_len = 50  # metres
    rct_tile_overlap = 10  # metres
    dir_matrix = Path('/data/TLS2/angola/2024-02-23_P02.RiSCAN/matrix/')
    fp_tile = Path('/data/TLS2/angola/2024-02-23_P02.RiSCAN/rct_tile50m_overlap10m/tile_index.dat')
    parent_dir = Path('/data/TLS2/angola/2024-02-23_P02.RiSCAN/rct_tile50m_overlap10m/tiled/')
    odir = Path('/data/TLS2/angola/results/tls_tree_attr/')
    figdir = Path('/data/TLS2/angola/results/figures/')
    num_workers = max(1, multiprocessing.cpu_count() - 1)
    
    # Validate output directories exist (for interactive mode)
    if not odir.exists():
        raise ValueError(f"Output directory does not exist: {odir}\nPlease create it first or specify an existing directory.")
    if not figdir.exists():
        raise ValueError(f"Figure directory does not exist: {figdir}\nPlease create it first or specify an existing directory.")
    
    # Setup logging for interactive mode
    logger = setup_logging(odir, country, plotid, date)

# ----------------------------------------------------------------------------
# Log configuration parameters
# ----------------------------------------------------------------------------
logger.info(f"Country: {country}")
logger.info(f"Plot ID: {plotid}")
logger.info(f"Date: {date}")
logger.info(f"RCT tile length: {rct_tile_len} m")
logger.info(f"RCT tile overlap: {rct_tile_overlap} m")
logger.info(f"Number of worker cores: {num_workers}")
logger.info("")

# ----------------------------------------------------------------------------
# Log input/output directories
# ----------------------------------------------------------------------------
logger.info("Input directories:")
logger.info(f"  SOP matrix dir: {dir_matrix}")
logger.info(f"  SOP matrix files found: {len(dat_files)}")
logger.info(f"  Tile index file: {fp_tile}")
# Log whether the tile index file exists
if fp_tile.exists():
    logger.info(f"  Tile index file exists: {fp_tile}")
else:
    logger.error(f"  Tile index file NOT found: {fp_tile}. Please provide a valid tile index file.")
    raise ValueError(f"Tile index file does not exist: {fp_tile}")
logger.info(f"  RCT extraction dir: {parent_dir}")
logger.info("")
logger.info("Output directories:")
logger.info(f"  CSV output dir: {odir}")
logger.info(f"  Figure output dir: {figdir}")
logger.info("")

# ----------------------------------------------------------------------------
# Discover input files
# ----------------------------------------------------------------------------
logger.info("Discovering input files...")
treeinfo_files = sorted(parent_dir.glob('*_trees_info.txt'))
tree_files = sorted(parent_dir.glob('*_trees.txt'))
leaves_files = sorted(parent_dir.glob('*/*_leaves.ply'))

# Convert Path objects to strings for compatibility with existing code
treeinfo_files = [str(f) for f in treeinfo_files]
tree_files = [str(f) for f in tree_files]
leaves_files = [str(f) for f in leaves_files]

logger.info(f"Found {len(treeinfo_files)} treeinfo files")
if treeinfo_files:
    logger.info(f"  Example: {treeinfo_files[0]}")
logger.info(f"Found {len(tree_files)} tree files")
if tree_files:
    logger.info(f"  Example: {tree_files[0]}")
logger.info(f"Found {len(leaves_files)} leaves files")
if leaves_files:
    logger.info(f"  Example: {leaves_files[0]}")
logger.info("")

#%%
# ============================================================================
# TREE DATA PARSING FUNCTIONS
# ============================================================================

# Attribution: Adapted from PyTreeFile (Tim Devereux)
# Source: https://github.com/tim-devereux/PyTreeFile/blob/main/pytreefile/treefiles.py
# Original author: Tim Devereux
# Notes: Function adapted to handle RayCloudTools `treeinfo` outputs and to
#        return a pandas.DataFrame consistent with this project's pipeline.

def treeinfo_attributes_tree(tree_file):
    """
    Adapted from: Tim Devereux / PyTreeFile
    Source: https://github.com/tim-devereux/PyTreeFile/blob/main/pytreefile/treefiles.py

    Extracts per-tree tree attributes from a tree file generated using treeinfo and returns a DataFrame.
    Can be run on a 'forest' or single treefile.

    Parameters:
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
    df = pd.DataFrame(line_list[1:], columns=line_list[0]).astype(float)
    if len(df) > 1:
        df.insert(0, "tree_id", range(1, len(df) + 1))

    return df


def treeinfo_attributes_segment(tree_file):
    """
    Adapted from: Tim Devereux / PyTreeFile
    Source: https://github.com/tim-devereux/PyTreeFile/blob/main/pytreefile/treefiles.py

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
    df = pd.DataFrame(line_list[1:], columns=line_list[0]).astype(float)
    df.insert(0, "tree_id", tree_ids)
    # remove row where parent_id is -1.0
    df = df[df["parent_id"] != -1.0]
    return df


def attributes_tree(tree_file):
    """
    Adapted from: Tim Devereux / PyTreeFile
    Source: https://github.com/tim-devereux/PyTreeFile/blob/main/pytreefile/treefiles.py

    Extracts per-tree tree attributes from a tree file generated using rayextract trees and returns a DataFrame.
    Can be run on a 'forest' or single treefile.

    Parameters:
    tree_file (str): The path to the treefile created using rayextract.

    Returns:
    pandas.DataFrame: A DataFrame containing tree attributes. If the input file is a forest, the DataFrame will contain a column 'tree_id' with the tree ID.
    """
    tree_ids = []
    line_list = []
    tree_id = 0

    with open(tree_file, "r") as file:
        lines = file.readlines()
        line_count = 0
        for line in lines[1:]:
            data = line.split(", ")
            for row in data:
                cell_data = data[0].strip().split(",")
                line_list.append(cell_data)
                tree_id += 1
            line_count += 1
    df = pd.DataFrame(line_list[2:], columns=line_list[0]).astype(float)
    # Remove duplicate rows
    df = df.drop_duplicates()

    if len(df) > 1:
        df.insert(0, "tree_id", range(1, len(df) + 1))

    return df

#%%
# ============================================================================
# STEP 1: READ AND MERGE TREE-LEVEL DATA
# ============================================================================
logger.info("STEP 1: Reading tree-level data from treeinfo files...")

## Read tree-level data from treeinfo.txt files
df_trees_info = pd.DataFrame()
for fn in treeinfo_files:
    # Extract tile coordinates from filename
    x, y = os.path.split(fn)[1].split('.')[0].split('_')[3:5]
    tile = f'Tile_{x}_{y}'
    
    # Parse tree attributes
    df = treeinfo_attributes_tree(fn)
    df['tile'] = tile
    
    # Move 'tile' column to the front
    cols = list(df.columns)
    cols = [cols[-1]] + cols[:-1]
    df = df[cols]
    df_trees_info = pd.concat([df_trees_info, df])

logger.info(f"Loaded {len(df_trees_info)} tree records from treeinfo files")
df_trees_info

#%%
# ============================================================================
# STEP 2: READ TREE LOCATION DATA
# ============================================================================
logger.info("STEP 2: Reading tree locations from trees.txt files...")

## Read tree locations from trees.txt files
df_trees = pd.DataFrame()
for fn in tree_files:
    # Extract tile coordinates from filename
    x, y = os.path.split(fn)[1].split('.')[0].split('_')[3:5]
    tile = f'Tile_{x}_{y}'
    
    # Parse tree attributes
    df = attributes_tree(fn)
    df['tile'] = tile
    
    # Move 'tile' column to the front
    cols = list(df.columns)
    cols = [cols[-1]] + cols[:-1]
    df = df[cols]
    df_trees = pd.concat([df_trees, df])

logger.info(f"Loaded {len(df_trees)} tree location records")
df_trees

#%%
# ============================================================================
# STEP 3: MERGE AND CLEAN TREE DATA
# ============================================================================
logger.info("STEP 3: Merging treeinfo and location data...")

## Merge treeinfo and trees dataframes
df_trees_merged = df_trees_info.merge(df_trees, on=['tile', 'tree_id'])
logger.info(f"Merged data: {len(df_trees_merged)} records")

# Drop duplicated trees in different tiles (keep unique trees only)
df_trees_merged = df_trees_merged.drop_duplicates(subset=['x', 'y', 'z'], keep=False)
logger.info(f"After removing duplicate trees across tiles: {len(df_trees_merged)} records")

# Add plot metadata
df_trees_merged['plot_id'] = plotid.split('_')[-1]

# Rearrange columns: put plot_id first, remove unnecessary columns
cols = list(df_trees_merged.columns)
cols = [cols[-1]] + cols[:-1]
df_trees_merged = df_trees_merged[cols]
df_trees_merged = df_trees_merged.drop(columns=['parent_id', 'section_id'])

logger.info(f"Final merged dataset: {len(df_trees_merged)} trees with {len(df_trees_merged.columns)} attributes")
df_trees_merged

#%%
# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

def plot_tiles(fp):
    """
    Visualize RCT tile layout with tile indices.
    
    Args:
        fp (Path or str): Path to tile index file (tile_index.dat)
        
    Returns:
        matplotlib.figure.Figure: Figure object with tile visualization
        
    Notes:
        - Each tile is drawn as a rectangle with its index labeled
        - Useful for understanding the tiling scheme
    """
    data = pd.read_csv(str(fp), sep=r'\s+', header=None, names=["Tile_Index", "Xmin", "Xmax", "Ymin", "Ymax"])
    fig, ax = plt.subplots(figsize=(8, 8))

    # Plot each tile as a square and add the tile index on top of each
    for _, row in data.iterrows():
        tile_index, xmin, xmax, ymin, ymax = row["Tile_Index"], row["Xmin"], row["Xmax"], row["Ymin"], row["Ymax"]
        # Plot the bounding box of the tile
        rect = plt.Rectangle((xmin, ymin), xmax - xmin, ymax - ymin, 
                             edgecolor='grey', facecolor='none', linewidth=1)
        ax.add_patch(rect)
 
        # Add tile index text in the centre of the tile
        ax.text((xmin + xmax) / 2., (ymin + ymax) / 2., str(tile_index), 
                ha='center', va='center', fontsize=8)

    # Set limits and aspect ratio
    ax.set_xlim(data["Xmin"].min() - 10, data["Xmax"].max() + 30)
    ax.set_ylim(data["Ymin"].min() - 10, data["Ymax"].max() + 30)
    ax.set_aspect('equal', adjustable='box')

    ax.set_xlabel("X Coordinate (m)", fontsize=12)
    ax.set_ylabel("Y Coordinate (m)", fontsize=12)
    ax.set_title(f"{date} {country} {plotid}", fontsize=12)

    plt.grid(visible=True, linestyle='--', alpha=0.5)
    plt.axis('equal')
    return fig


# fig = plot_tiles(fp_tile)

# %%
def plot_tiles_scans_boundary(fp_tile, dir_matrix, show=True,
                              savefig=False, figfn=None):
    """
    Create comprehensive visualization of scan positions, plot boundary, and tile layout.
    
    This function generates a plot showing:
    - Scan positions (blue circles)
    - Plot boundary (red dashed line - minimum rotated rectangle)
    - Tile grid (gray rectangles with indices)
    
    Args:
        fp_tile (Path or str): Path to tile index file
        dir_matrix (Path or str): Path to directory containing SOP matrix files
        show (bool): Whether to display the plot (default: True)
        savefig (bool): Whether to save the figure to file (default: False)
        figfn (Path or str): Output filename if savefig=True
        
    Returns:
        matplotlib.axes.Axes: Axes object with the plot
        
    Raises:
        ValueError: If no .DAT files found in dir_matrix
        ValueError: If savefig=True but figfn is None
        
    Notes:
        - Plot area is calculated as minimum rotated rectangle around scan positions
        - Scan positions are labeled with their ScanPos ID
    """
    # Convert to Path objects for robust path handling
    fp_tile = Path(fp_tile)
    dir_matrix = Path(dir_matrix)
    
    # Read all scan positions
    scan_positions = []
    scan_ids = []

    # Look for both .DAT and .dat files (case-insensitive)
    dat_files = list(dir_matrix.glob('*.DAT')) + list(dir_matrix.glob('*.dat'))
    if not dat_files:
        raise ValueError(f"No .DAT or .dat files found in {dir_matrix}. Please check the directory path.")

    for dat in dat_files:
        scan_id = dat.stem.split('ScanPos')[1]
        scan_coord = np.loadtxt(str(dat))[:3, 3]
        scan_positions.append(scan_coord)
        scan_ids.append(scan_id)

    # Convert to DataFrame and then to GeoDataFrame
    df_sp = pd.DataFrame(scan_positions, columns=['x', 'y', 'z'])
    df_sp['ScanPos'] = scan_ids
    df_sp['geometry'] = [Point(x, y) for x, y in zip(df_sp['x'], df_sp['y'])]

    # Create GeoDataFrame and explicitly set geometry
    sp = gp.GeoDataFrame(df_sp, geometry='geometry')

    # Compute and print the plot area
    area = sp.unary_union.minimum_rotated_rectangle.area
    logger.info(f'Plot area: {area / 1e4:.2f} ha')

    # Plot scan locations
    ax = sp.plot(figsize=(8, 8), marker='o', markersize=20,
                 color='b', alpha=0.3, label='scan locations')

    # Add the ScanPos as text above each point
    for j, row in sp.iterrows():
        spn = int(row.ScanPos)
        ha = 'right' if spn % 2 == 0 else 'left'
        va = 'top' if spn % 2 == 0 else 'bottom'
        ax.text(row.x, row.y, row.ScanPos, fontsize=8, 
                ha=ha, va=va, color='b', alpha=0.8)

    # Plot boundary
    boundary = sp.unary_union.minimum_rotated_rectangle
    ax.plot(*boundary.exterior.xy, color='r', lw=1, ls='--', 
            zorder=1, label='plot boundary')

    # Read tile index
    data = pd.read_csv(str(fp_tile), sep=r'\s+', header=None, 
                       names=["Tile_Index", "Xmin", "Xmax", "Ymin", "Ymax"])

    # Plot each tile as a square and add the tile index on top of each
    for _, row in data.iterrows():
        tile_index, xmin, xmax, ymin, ymax = row["Tile_Index"], row["Xmin"], row["Xmax"], row["Ymin"], row["Ymax"]
        # Plot the bounding box of the tile
        rect = plt.Rectangle((xmin, ymin), xmax - xmin, ymax - ymin, 
                             edgecolor='grey', facecolor='none', linewidth=1)
        ax.add_patch(rect)
 
        # Add tile index text in the centre of the tile
        ax.text((xmin + xmax) / 2., (ymin + ymax) / 2., str(tile_index), 
                ha='center', va='center', fontsize=8)

    # Set limits and aspect ratio
    ax.set_xlim(data["Xmin"].min() - 10, data["Xmax"].max() + 30)
    ax.set_ylim(data["Ymin"].min() - 10, data["Ymax"].max() + 30)
    ax.set_aspect('equal', adjustable='box')

    ax.set_xlabel("X Coordinate (m)", fontsize=12)
    ax.set_ylabel("Y Coordinate (m)", fontsize=12)

    ax.set_title(f"{date} {country} {plotid} (plot area: {area / 1e4:.2f} ha)", 
                 fontsize=12)
    ax.legend(loc='upper right', fontsize=10)
    ax.axis('equal')

    plt.grid(visible=True, linestyle='--', alpha=0.5)
    if savefig:
        if figfn is None:
            raise ValueError("figfn must be provided when savefig=True")
        figfn = Path(figfn)
        figfn.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(str(figfn), dpi=300)
        logger.info(f'Figure saved to {figfn}')

    if show:
        plt.show()
    else:
        plt.close()
    
    return ax

### temp test
# fig = plot_tiles_scans_boundary(fp_tile, dir_matrix, show=True)

# %%
# ============================================================================
# STEP 4: GENERATE PLOT VISUALIZATION
# ============================================================================
logger.info("STEP 4: Generating plot visualization with scan positions and tiles...")

### save figure
tmp_ofn = figdir / f'{date}_{country}_{plotid}_scanloc_rct_tile{rct_tile_len}m_overlap{rct_tile_overlap}m.png'
plot_tiles_scans_boundary(fp_tile, dir_matrix, show=False, savefig=True, figfn=tmp_ofn)
logger.info("")

#%%
# ============================================================================
# STEP 5: SPATIAL FILTERING - IDENTIFY IN-PLOT TREES
# ============================================================================
logger.info("STEP 5: Calculating plot boundary and filtering in-plot trees...")

df = df_trees_merged.copy()

# Calculate plot boundary based on scan positions
sp = gp.GeoDataFrame(columns=['x', 'y', 'z', 'sp'], geometry=[])

# Use Path for robust directory handling
dir_matrix_path = Path(dir_matrix)
# Look for both .DAT and .dat files (case-insensitive)
dat_files = list(dir_matrix_path.glob('*.DAT')) + list(dir_matrix_path.glob('*.dat'))
if not dat_files:
    raise ValueError(f"No .DAT or .dat files found in {dir_matrix_path}. Please check the directory path.")

for i, dat in enumerate(dat_files):
    sp.loc[i, 'ScanPos'] = dat.stem.split('ScanPos')[1]
    sp.iloc[i, :3] = np.loadtxt(str(dat))[:3, 3]

# Set geometry explicitly to avoid warnings
sp['geometry'] = [Point(r.x, r.y) for r in sp.itertuples()]
sp = sp.set_geometry('geometry')

# Compute plot boundary (minimum rotated rectangle around scan positions)
extent = sp.unary_union.minimum_rotated_rectangle
area = extent.area
logger.info(f'Plot area (minimum rotated rectangle): {area / 1e4:.2f} ha')

# Find trees within the plot boundary
df['geometry'] = [Point(x, y) for x, y in zip(df['x'], df['y'])]
df = gp.GeoDataFrame(df, geometry='geometry')  # Explicitly setting the geometry
df['in_plot'] = df.within(extent)
logger.info(f'Total trees extracted: {len(df)}')
logger.info(f'In-plot trees: {df["in_plot"].sum()}')
logger.info(f'Out-of-plot trees: {(~df["in_plot"]).sum()}')

# Add metadata columns
df['date'] = date
df['country'] = country

# Rearrange columns: date, country, plot_id, tile, tree_id, in_plot, ...
cols = list(df.columns)
cols = [cols[-2]] + [cols[-1]] + cols[:3] + ['in_plot'] + cols[3:-3]
df = df[cols]
df

#%%
## Filter trees of interest (in-plot, with minimum height and DBH thresholds)
trees_interested = df[(df.in_plot == True) & (df.height >= 3) & (df.DBH >= 0.1)]
logger.info(f'RCT extraction results: {len(trees_interested)} in-plot trees with H >= 3m and DBH >= 10cm')

# %%
# ============================================================================
# LEAF AREA CALCULATION FUNCTIONS
# ============================================================================

def leaf_stats(leaf_file):
    """
    Calculate leaf area from triangulated leaf mesh files.
    
    This function reads PLY mesh files containing triangulated leaf surfaces
    and calculates the total one-sided leaf area by summing triangle areas.
    
    Args:
        leaf_file (str): Path to PLY file containing leaf mesh (*_leaves.ply)
        
    Returns:
        tuple: (num_leaves, total_area, tri_areas)
            - num_leaves (int): Number of triangles (leaf count)
            - total_area (float): Total one-sided leaf area (m²)
            - tri_areas (np.ndarray): Individual triangle areas (m²)
            
    Notes:
        - Leaf area is calculated using cross product of triangle edge vectors
        - Formula: Area = 0.5 * ||(v1-v0) × (v2-v0)||
        - This gives one-sided leaf area (not total leaf surface area)
    """
    mesh = PlyData.read(leaf_file)
    # Extract coordinates of leaf points (vertices)
    vertices = np.vstack([mesh['vertex'][axis] for axis in ['x', 'y', 'z']]).T
    # Extract vertex indices to construct triangle polygons for each leaf
    faces = np.vstack(mesh['face']['vertex_indices']) if 'face' in mesh else None
    v0 = vertices[faces[:,0]]
    v1 = vertices[faces[:,1]]
    v2 = vertices[faces[:,2]]
    # Count the number of leaves (triangles)
    num_leaves = faces.shape[0]
    # Calculate area of each triangle: 0.5 * || (v1 - v0) × (v2 - v0) ||
    tri_areas = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
    # Sum to get total one-side leaf area
    total_area = tri_areas.sum()

    return num_leaves, total_area, tri_areas

# ### temp test
# leaf_fn = leaves_files[0]
# plot_id = os.path.basename(leaf_fn).split('_')[1]
# tx, ty = os.path.basename(leaf_fn).split('_')[3:5]
# tile_id = f'Tile_{tx}_{ty}'
# tree_id = os.path.basename(leaf_fn).split('_')[-2]

# num_leaves, total_area, tri_areas = leaf_stats(leaf_fn)
# print(f'Plot: {plot_id}, {tile_id}, Tree_id: {tree_id}')
# print(f'Number of leaves (triangles): {num_leaves}')
# print(f'Total one-side leaf area: {total_area:.2f} m²')

#%%
## parallelized reading leaf mesh files
def process_leaf_file(leaf_fn):
    """
    Worker function for parallel processing of leaf mesh files.
    
    Extracts leaf statistics from a single PLY file. This function is designed
    to be called in parallel by ProcessPoolExecutor.
    
    Args:
        leaf_fn (str): Path to leaf mesh PLY file
        
    Returns:
        dict: Dictionary containing:
            - plot_id: Plot ID extracted from filename
            - tile_id: Tile ID extracted from filename
            - tree_id: Tree ID extracted from filename  
            - num_leaves: Number of leaf triangles
            - leaf_area_m2: Total leaf area in m²
            - error: Error message if processing failed (None if successful)
            
    Notes:
        - Defined at module level so it can be pickled by ProcessPoolExecutor
        - Returns error dict with zeros if processing fails
    """
    try:
        # Extract metadata from filename
        base = os.path.basename(leaf_fn).split('_')
        plot_id = base[1] if len(base) > 1 else None
        tx, ty = (base[3], base[4]) if len(base) > 4 else (None, None)
        tile_id = f'Tile_{tx}_{ty}' if tx is not None and ty is not None else None
        tree_id = base[-2] if len(base) >= 2 else None

        # Calculate leaf statistics
        num_leaves, total_area, _ = leaf_stats(leaf_fn)

        return {
            'plot_id': plot_id,
            'tile_id': tile_id,
            'tree_id': tree_id,
            'num_leaves': int(num_leaves),
            'leaf_area_m2': round(float(total_area), 3),
            'error': None
        }
    except Exception as e:
        return {
            'plot_id': None,
            'tile_id': None,
            'tree_id': None,
            'num_leaves': 0,
            'leaf_area_m2': 0.0,
            'error': str(e)
        }

#%%
# ============================================================================
# STEP 6: PARALLEL PROCESSING OF LEAF MESH FILES
# ============================================================================
logger.info("STEP 6: Processing leaf mesh files to calculate leaf areas...")
logger.info(f"Using {num_workers} parallel workers")

# Process all leaf files in parallel
rows = []

with ProcessPoolExecutor(max_workers=num_workers) as exe:
    futures = {exe.submit(process_leaf_file, fn): fn for fn in leaves_files}
    # Progress bar to track completion
    for fut in tqdm(as_completed(futures), total=len(futures), desc="Processing leaf files"):
        res = fut.result()
        rows.append(res)

# Build DataFrame from results
df_leaf = pd.DataFrame(rows)

# Report processing results
n_errors = df_leaf['error'].notnull().sum()
logger.info(f'Completed processing {len(df_leaf)} leaf files')
logger.info(f'Successful: {len(df_leaf) - n_errors}, Errors: {n_errors}')

# Keep only useful columns and convert tree_id to int to match main df
df_leaf = df_leaf[['plot_id','tile_id','tree_id','num_leaves','leaf_area_m2','error']]
df_leaf['tree_id'] = df_leaf['tree_id'].astype(int)

df_leaf

# %%
# ============================================================================
# STEP 7: MERGE LEAF DATA WITH TREE ATTRIBUTES
# ============================================================================
logger.info("STEP 7: Merging leaf area data with tree attributes...")

# Append leaf data columns to main df by merging on tile and tree_id
df = df.merge(df_leaf[['tile_id', 'tree_id', 'num_leaves', 'leaf_area_m2']], 
              left_on=['tile', 'tree_id'], 
              right_on=['tile_id', 'tree_id'], 
              how='left')
df = df.drop(columns=['tile_id'])

# Fill NaN values with 0 (trees without leaf data)
df['num_leaves'] = df['num_leaves'].fillna(0).astype(int)
df['leaf_area_m2'] = df['leaf_area_m2'].fillna(0).astype(float)

trees_with_leaves = (df['leaf_area_m2'] > 0).sum()
logger.info(f'Trees with leaf area data: {trees_with_leaves} / {len(df)}')
logger.info("")

df

#%%
# ============================================================================
# STEP 8: EXPORT RESULTS TO CSV
# ============================================================================
logger.info("STEP 8: Saving results to CSV...")

# Generate output filename
ofn = odir / f'{date}_{country}_{plotid}_treeLevel_attr_rct_tile{rct_tile_len}m_overlap{rct_tile_overlap}m.csv'
logger.info(f'Output file: {ofn}')

# Save to CSV
df.to_csv(str(ofn), index=False)
logger.info(f'Saved {len(df)} tree records with {len(df.columns)} attributes')

# Log summary statistics
logger.info("")
logger.info("="*80)
logger.info("PROCESSING COMPLETE - SUMMARY STATISTICS")
logger.info("="*80)
logger.info(f"Total trees extracted: {len(df)}")
logger.info(f"In-plot trees: {df['in_plot'].sum()}")
logger.info(f"Trees with H>=3m & DBH>=10cm: {len(trees_interested)}")
logger.info(f"Trees with leaf data: {trees_with_leaves}")
logger.info(f"")
logger.info(f"Output CSV: {ofn}")
logger.info(f"Output figure: {tmp_ofn}")
logger.info("="*80)
logger.info("Pipeline completed successfully!")
logger.info("="*80)
