import matplotlib
matplotlib.use("Agg")  # Use non-GUI backend for image generation
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, MultiPolygon
from matplotlib.patches import Polygon as mplPolygon
from PIL import Image
import math
import json
import numpy as np
import time
import openslide
from openslide import OpenSlide
from collections import defaultdict

# Note on GPU acceleration:
# While the current implementation uses Shapely (CPU-bound),
# libraries like cuSpatial could be used for GPU acceleration
# of spatial operations if available in your environment.
# This would require a more substantial rewrite.

# Define QC class mapping (for reference, not directly used by drawing if names are in JSON)
QC_CLASS_MAPPING = {
    1: "Normal Tissue",
    2: "Fold",
    3: "Darkspot & Foreign Object",
    4: "PenMarking",
    5: "Edge & Air Bubble",
    6: "OOF",  # Out of Focus
    7: "Background"
}

# Define colors for QC annotation types
# TLS masks are 'lime' (approx #00FF00)
QC_COLOR_MAPPING = {
    "Fold":                         '#E6194B',  # Red
    "PenMarking":                   '#4363D8',  # Blue
    "Edge & Air Bubble":            '#F58231',  # Orange
    "OOF":                          '#911EB4',  # Purple
    "Darkspot & Foreign Object":    '#FFE119',  # Yellow
    "Normal Tissue":                '#A9A9A9',  # DarkGray (for less emphasis)
    "Background":                   '#D3D3D3',  # LightGray (for even less emphasis)
    # Default color if type not found, can be handled in the drawing function
}

def _precompute_tls_qc_intersections(tls_polygons, qc_polygons_with_type):
    """
    Pre-compute all intersections between TLS and QC polygons.
    
    Args:
        tls_polygons (list): List of Shapely polygons for TLS annotations.
        qc_polygons_with_type (list): List of dicts with QC polygons and their types.
        
    Returns:
        dict: A dictionary mapping TLS polygon indices to lists of intersecting QC polygons.
    """
    intersection_map = defaultdict(list)
    
    for tls_idx, tls_poly in enumerate(tls_polygons):
        if tls_poly is None or tls_poly.is_empty:
            continue
            
        for qc_info in qc_polygons_with_type:
            qc_poly = qc_info['polygon']
            
            if qc_poly is None or qc_poly.is_empty:
                continue
                
            if tls_poly.intersects(qc_poly):
                intersection_map[tls_idx].append(qc_info)
    
    return intersection_map

def _draw_overlapping_qc_annotations(ax, intersecting_qc_polygons, patch_top_left_coords_level0, qc_color_map, downsample_factor=1):
    """
    Draws QC annotations that overlap with the current TLS polygon onto the given axes.
    
    Args:
        ax: Matplotlib axes object to draw on.
        intersecting_qc_polygons (list): List of QC polygons known to intersect with the current TLS.
        patch_top_left_coords_level0 (tuple): (x, y) coordinates of the top-left corner of the patch.
        qc_color_map (dict): Mapping from QC type name to color string.
        downsample_factor (float): Downsample factor for scaling QC polygons.
    """
    for qc_info in intersecting_qc_polygons:
        qc_shapely_poly = qc_info['polygon']  # This is at level 0
        qc_type_name = qc_info['type_name']
        
        qc_color = qc_color_map.get(qc_type_name, '#808080')  # Default color
        
        geoms_to_draw = []
        if isinstance(qc_shapely_poly, Polygon):
            geoms_to_draw.append(qc_shapely_poly)
        elif isinstance(qc_shapely_poly, MultiPolygon):
            geoms_to_draw.extend(list(qc_shapely_poly.geoms))
        
        for geom in geoms_to_draw:
            # Translate coordinates relative to the patch's top-left (at level 0), then scale to current level
            translated_scaled_coords = [
                ((x - patch_top_left_coords_level0[0]) / downsample_factor,
                 (y - patch_top_left_coords_level0[1]) / downsample_factor)
                for x, y in geom.exterior.coords
            ]
            
            mpl_qc_poly = mplPolygon(translated_scaled_coords, fill=True, color=qc_color, alpha=0.4, linewidth=1, edgecolor=qc_color)
            ax.add_patch(mpl_qc_poly)

def tls_annotations_json(json_data):
    """
    Load annotations from JSON Object.
    
    Args:
        json_data (str): JSON Data object
        
    Returns:
        dict: Dictionary containing annotation data.
    """
    annotations = json_data

    num_annotations = len(annotations['features'])
    all_annotations = {}


    if num_annotations == 0:
        raise ValueError("No annotations found in the JSON file.")

    for feature_num in range(num_annotations):
        name = annotations['features'][feature_num]['properties']['name'] + ": TLS number " + str(feature_num)
        all_annotations[name] = {}
        all_annotations[name]['coordinates'] = annotations['features'][feature_num]['geometry']['coordinates'][0]

    return all_annotations


def qc_annotations_json(json_data):
    """
    Load annotations from JSON Object.
    
    Args:
        json_data (str): JSON Data object
        
    Returns:
        dict: Dictionary containing annotation data.
    """
    annotations = json_data

    num_annotations = len(annotations['features'])
    all_annotations = {}


    if num_annotations == 0:
        raise ValueError("No annotations found in the JSON file.")

    for feature_num in range(num_annotations):
        name =  "# " + str(feature_num)

        type = annotations['features'][feature_num]['properties']['classification']

        all_annotations[name] = {}
        all_annotations[name]['coordinates'] = annotations['features'][feature_num]['geometry']['coordinates'][0]
        all_annotations[name]['type'] = type
        all_annotations[name]['area'] = annotations['features'][feature_num]['properties']['area']

    return all_annotations

def create_polygon_from_coordinates(coordinates):
    """
    Create a shapely Polygon from coordinates.
    
    Args:
        coordinates (list): List of [x, y] coordinates.
        
    Returns:
        Polygon: A shapely Polygon object.
    """
    return Polygon(coordinates)

def _load_and_prepare_annotations(tls_annotations, qc_annotations):
    """
    Load and prepare TLS and QC annotations from dict data.
    
    Args:
        tls_annotations (dict): Dictionary containing TLS annotation data.
        qc_annotations (dict): Dictionary containing QC annotation data.
        
    Returns:
        tuple: (all_tls_polygons, all_tls_labels, all_qc_polygons_with_type)
            - all_tls_polygons: List of Shapely polygons for TLS annotations
            - all_tls_labels: List of labels for TLS annotations
            - all_qc_polygons_with_type: List of dicts with QC polygons and their types
    """
    all_tls_polygons = []
    all_tls_labels = []
    
    # Process TLS annotations
    for annotation_id, data in tls_annotations.items():
        polygon = None
        if 'coordinates' in data and data['coordinates']:
            try:
                polygon = create_polygon_from_coordinates(data['coordinates'])
            except Exception as e:
                print(f"Warning: Could not create polygon for {annotation_id}: {e}")
        
        all_tls_labels.append(annotation_id)
        all_tls_polygons.append(polygon)
    
    # Process QC annotations
    all_qc_polygons_with_type = []
    if qc_annotations:
        for qc_id, qc_data in qc_annotations.items():
            qc_poly = None
            if 'coordinates' in qc_data and qc_data['coordinates']:
                try:
                    qc_poly = create_polygon_from_coordinates(qc_data['coordinates'])
                except Exception as e:
                    print(f"Warning: Could not create QC polygon for {qc_id}: {e}")
            if qc_poly and not qc_poly.is_empty:
                all_qc_polygons_with_type.append({
                    'polygon': qc_poly,
                    'type_name': qc_data['type']
                })
    
    return all_tls_polygons, all_tls_labels, all_qc_polygons_with_type

def _calculate_patch_parameters(tls_polygon, intersecting_qc_polygons, padding_factor=0.1):
    """
    Calculate parameters for extracting a patch from the WSI.
    
    Args:
        tls_polygon: Shapely polygon for the TLS annotation.
        intersecting_qc_polygons: List of QC polygons that intersect with this TLS.
        padding_factor: Factor to determine padding amount (as a fraction of the side length).
        
    Returns:
        tuple: (final_patch_location, patch_size_at_level0, side_length)
            - final_patch_location: Tuple (x, y) for top-left corner of patch (at level 0)
            - patch_size_at_level0: Size of the patch at level 0 
            - side_length: Side length of the square patch (at level 0)
    """
    # Get bounds of the TLS polygon (level 0)
    tls_minx, tls_miny, tls_maxx, tls_maxy = tls_polygon.bounds
    
    # Initialize combined bounds with the TLS polygon's bounds
    combined_minx, combined_miny, combined_maxx, combined_maxy = tls_minx, tls_miny, tls_maxx, tls_maxy
    
    # Expand bounds to include intersecting QC annotations
    for qc_info in intersecting_qc_polygons:
        qc_poly_level0 = qc_info['polygon']
        if not qc_poly_level0.is_empty:
            qc_bounds = qc_poly_level0.bounds
            combined_minx = min(combined_minx, qc_bounds[0])
            combined_miny = min(combined_miny, qc_bounds[1])
            combined_maxx = max(combined_maxx, qc_bounds[2])
            combined_maxy = max(combined_maxy, qc_bounds[3])
    
    # Calculate dimensions based on the combined bounds (at level 0)
    orig_patch_width = int(combined_maxx - combined_minx)
    orig_patch_height = int(combined_maxy - combined_miny)
    
    if orig_patch_width <= 0 or orig_patch_height <= 0:
        return None, None, None
    
    # Determine the side length for the square patch (at level 0)
    side_length = max(orig_patch_width, orig_patch_height)
    
    # Calculate the center of the combined bounding box (at level 0)
    center_x = combined_minx + orig_patch_width / 2.0
    center_y = combined_miny + orig_patch_height / 2.0
    
    # Determine the top-left corner of the square region (before padding, at level 0)
    square_region_minx = center_x - side_length / 2.0
    square_region_miny = center_y - side_length / 2.0
    
    # Add padding to this square region (padding amount at level 0)
    padding_amount = int(side_length * padding_factor)
    
    # Location for wsi.read_region (at level 0)
    read_location_x = int(square_region_minx - padding_amount)
    read_location_y = int(square_region_miny - padding_amount)
    # Total dimension of the square region to read, at level 0
    read_size_dimension = int(side_length + 2 * padding_amount)
    
    # Ensure location is not negative (level 0 coordinates for read_region)
    final_patch_location = (max(0, read_location_x), max(0, read_location_y))
    
    return final_patch_location, read_size_dimension, side_length

def _draw_single_patch_content(ax, tls_polygon, intersecting_qc_polygons, patch_np, final_patch_location, 
                              downsample_factor, label_text, with_mask=True):
    """
    Draw the content (TLS polygons, QC polygons, label) on a patch image.
    
    Args:
        ax: Matplotlib axes to draw on.
        tls_polygon: Shapely polygon for the TLS annotation.
        intersecting_qc_polygons: List of QC polygons that intersect with this TLS.
        patch_np: NumPy array representation of the patch image.
        final_patch_location: Tuple (x, y) for top-left corner of patch (at level 0).
        downsample_factor: Downsample factor for the current WSI level.
        label_text: Text label to display on the patch.
        with_mask: Whether to display polygon masks.
    """
    ax.imshow(patch_np)
    
    # Draw TLS mask first (if with_mask is True)
    if with_mask and tls_polygon:
        if isinstance(tls_polygon, Polygon):
            # Translate polygon coordinates to be relative to the final_patch_location (level 0), then scale to current level
            translated_coords = [
                ((x - final_patch_location[0]) / downsample_factor, 
                 (y - final_patch_location[1]) / downsample_factor) 
                for x, y in tls_polygon.exterior.coords
            ]
            mpl_tls_poly = mplPolygon(translated_coords, fill=False, color='lime', linewidth=2)
            ax.add_patch(mpl_tls_poly)
        elif isinstance(tls_polygon, MultiPolygon):
            for poly_geom in tls_polygon.geoms:
                translated_coords = [
                    ((x - final_patch_location[0]) / downsample_factor, 
                     (y - final_patch_location[1]) / downsample_factor) 
                    for x, y in poly_geom.exterior.coords
                ]
                mpl_tls_poly = mplPolygon(translated_coords, fill=False, color='lime', linewidth=2)
                ax.add_patch(mpl_tls_poly)
    
    # Draw overlapping QC annotations
    if tls_polygon and intersecting_qc_polygons:
        _draw_overlapping_qc_annotations(
            ax,
            intersecting_qc_polygons,
            final_patch_location, # Top-left of the patch at level 0
            QC_COLOR_MAPPING,
            downsample_factor 
        )
    
    # Add annotation label
    x_pos = patch_np.shape[1] * 0.05
    y_pos = patch_np.shape[0] * 0.1
    ax.text(x_pos, y_pos, label_text, fontsize=12, color='white',
            bbox=dict(facecolor='black', alpha=0.5, boxstyle='round,pad=0.3'))

def display_annotations(qc_annotations, tls_annotations, image_path, save_path=None, with_mask=True, 
                       num_per_row=7, file_format="png", level=0):
    """
    Visualize annotations with optional masks by extracting patches from WSI.
    
    Args:
        qc_annotations (dict): Dict containing qc_annotations (keys: 'type', 'coordinates').
        tls_annotations (dict): Dict containing tls_annotations and path coordinates.
        image_path (str): Path to the WSI image.
        save_path (str): Path to save the figure. If None, shows the plot.
        with_mask (bool): Whether to overlay mask polygons.
        num_per_row (int): Number of tls_annotations per row.
        file_format (str): Format for saving the figure.
        level (int): The resolution level to read from the WSI (0 is highest).
    """
    start_time = time.time()
    
    try:
        wsi = OpenSlide(image_path)
    except Exception as e:
        print(f"Error opening WSI {image_path}: {e}")
        # Create a dummy plot indicating the error
        fig, ax = plt.subplots(figsize=(4, 4))
        ax.text(0.5, 0.5, f"Error loading WSI:\n{image_path}", fontsize=10, ha="center", va="center", color="red")
        ax.axis('off')
        if save_path:
            plt.savefig(save_path, format=file_format, bbox_inches='tight')
        else:
            plt.show()
        plt.close()
        return

    if level < 0 or level >= wsi.level_count:
        raise ValueError(f"Invalid level {level}. Must be between 0 and {wsi.level_count - 1}.")

    downsample_factor = wsi.level_downsamples[level]
    print(f"Using WSI level {level} with downsample factor {downsample_factor}")
    
    # Load and prepare annotations (TLS and QC)
    all_tls_polygons, all_tls_labels, all_qc_polygons_with_type = _load_and_prepare_annotations(tls_annotations, qc_annotations)
    print(f"Loaded {len(all_tls_polygons)} TLS annotations and {len(all_qc_polygons_with_type)} QC annotations")
    
    # Pre-compute TLS-QC intersections
    tls_qc_intersection_map = _precompute_tls_qc_intersections(all_tls_polygons, all_qc_polygons_with_type)
    print(f"Pre-computed intersections in {time.time() - start_time:.2f} seconds")
    
    num_tls = len(all_tls_polygons)
    
    if num_tls == 0:
        fig, ax = plt.subplots(figsize=(4, 4))
        ax.text(0.5, 0.5, "No TLS present", fontsize=16, ha="center", va="center")
        ax.axis('off')
        if save_path:
            plt.savefig(save_path, format=file_format, bbox_inches='tight')
        else:
            plt.show()
        plt.close()
        return
    
    num_rows = math.ceil(num_tls / num_per_row)
    fig, axes = plt.subplots(num_rows, num_per_row, figsize=(num_per_row * 4, num_rows * 4))
    axes = axes.flatten() if num_tls > 1 else [axes]
    
    # Process each TLS annotation
    for i in range(num_tls):
        current_polygon = all_tls_polygons[i]
        ax = axes[i]
        ax.axis('off')  # Turn off axis for each subplot initially

        if current_polygon is None or current_polygon.is_empty:
            ax.text(0.5, 0.5, "No polygon data", fontsize=10, ha="center", va="center", transform=ax.transAxes)
            continue

        try:
            # Get pre-computed intersecting QC polygons for this TLS
            intersecting_qc_polygons = tls_qc_intersection_map.get(i, [])
            
            # Calculate patch parameters (location, size, etc.)
            final_patch_location, read_size_dimension, side_length = _calculate_patch_parameters(
                current_polygon, intersecting_qc_polygons
            )
            
            if final_patch_location is None or read_size_dimension is None or side_length is None:
                ax.text(0.5, 0.5, "Invalid combined area", fontsize=10, ha="center", va="center", transform=ax.transAxes)
                continue

            # Calculate the size of the region to read at the specified WSI level
            patch_read_size_at_level = (
                max(1, int(read_size_dimension / downsample_factor)),
                max(1, int(read_size_dimension / downsample_factor))
            )

            # Extract patch from WSI
            patch_pil = wsi.read_region(final_patch_location, level, patch_read_size_at_level)
            patch_np = np.array(patch_pil.convert("RGB"))
            
            # Prepare label text
            annotation_text = all_tls_labels[i]
            if len(annotation_text) > 20:
                mid = len(annotation_text) // 2
                # Ensure split doesn't happen at the very beginning or end if mid is 0 or len
                if mid > 0 and mid < len(annotation_text):
                    label_text = annotation_text[:mid] + "\n" + annotation_text[mid:]
                else:
                    label_text = annotation_text
            else: 
                label_text = annotation_text
            
            # Draw patch content (TLS mask, QC overlays, label)
            _draw_single_patch_content(
                ax, current_polygon, intersecting_qc_polygons, patch_np, 
                final_patch_location, downsample_factor, label_text, with_mask
            )
            
        except Exception as e:
            print(f"Error processing annotation {all_tls_labels[i]}: {e}")
            ax.text(0.5, 0.5, "Error processing", fontsize=10, color="red", ha="center", va="center", transform=ax.transAxes)

    # Turn off any extra axes
    for j in range(num_tls, len(axes)):
        axes[j].axis('off')
        
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, format=file_format, bbox_inches='tight')
    else:
        plt.show()
    plt.close()
    wsi.close()  # Close the WSI object
    
    end_time = time.time()
    print(f"Total processing time: {end_time - start_time:.2f} seconds")
