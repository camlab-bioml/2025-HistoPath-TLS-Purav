import json
import numpy as np
import time
import sys
import os
import csv
from shapely.geometry import Polygon, MultiPolygon
from collections import defaultdict
import openslide
from openslide import OpenSlide

# Define QC class mapping
QC_CLASS_MAPPING = {
    2: "Fold",
    3: "Darkspot & Foreign Object", 
    4: "PenMarking",
    5: "Edge & Air Bubble",
    6: "OOF"  # Out of Focus
}

def extract_mpp_from_wsi(wsi_path):
    """
    Extract MPP (Microns Per Pixel) values from WSI metadata.
    
    Args:
        wsi_path (str): Path to the WSI file
        
    Returns:
        tuple: (mpp_x, mpp_y) or (None, None) if not available
    """
    try:
        wsi = OpenSlide(wsi_path)
        properties = wsi.properties
        
        # Try different property names used by different vendors
        mpp_x = None
        mpp_y = None
        
        # Aperio/Leica format
        if 'aperio.MPP' in properties:
            mpp_x = mpp_y = float(properties['aperio.MPP'])
        elif 'openslide.mpp-x' in properties and 'openslide.mpp-y' in properties:
            mpp_x = float(properties['openslide.mpp-x'])
            mpp_y = float(properties['openslide.mpp-y'])
        # Try alternative property names
        elif 'tiff.XResolution' in properties and 'tiff.YResolution' in properties:
            # Convert from resolution to MPP (assuming resolution is in pixels per cm)
            try:
                x_res = float(properties['tiff.XResolution'])
                y_res = float(properties['tiff.YResolution'])
                if x_res > 0 and y_res > 0:
                    mpp_x = 10000.0 / x_res  # Convert pixels/cm to microns/pixel
                    mpp_y = 10000.0 / y_res
            except (ValueError, ZeroDivisionError):
                pass
        
        wsi.close()
        
        if mpp_x is not None and mpp_y is not None:
            print(f"Extracted MPP from {wsi_path}: X={mpp_x:.4f}, Y={mpp_y:.4f} microns/pixel")
            return mpp_x, mpp_y
        else:
            print(f"Warning: Could not extract MPP from {wsi_path}")
            return None, None
            
    except Exception as e:
        print(f"Error reading WSI {wsi_path}: {e}")
        return None, None

def convert_area_pixels_to_microns(area_pixels, mpp_x, mpp_y):
    """
    Convert area from pixels² to microns².
    
    Args:
        area_pixels (float): Area in pixels²
        mpp_x (float): Microns per pixel in X direction
        mpp_y (float): Microns per pixel in Y direction
        
    Returns:
        float: Area in microns²
    """
    if mpp_x is None or mpp_y is None:
        return None
    return area_pixels * mpp_x * mpp_y

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


def _create_polygon_from_coordinates(coordinates):
    """
    Create a shapely Polygon from coordinates with validation.
    
    Args:
        coordinates (list): List of [x, y] coordinates.
        
    Returns:
        Polygon: A shapely Polygon object, or None if invalid.
    """
    try:
        polygon = Polygon(coordinates)
        
        # Check if polygon is valid
        if not polygon.is_valid:
            # Try to fix the polygon using buffer(0) - common fix for self-intersections
            fixed_polygon = polygon.buffer(0)
            if fixed_polygon.is_valid and fixed_polygon.area > 0:
                return fixed_polygon
            else:
                # If we can't fix it, return None
                return None
        
        # Check for reasonable area (not zero or negative)
        if polygon.area <= 0:
            return None
            
        return polygon
    except Exception as e:
        # Silently return None for invalid coordinates
        return None

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
    tls_invalid_count = 0
    
    # Process TLS annotations
    for annotation_id, data in tls_annotations.items():
        polygon = None
        if 'coordinates' in data and data['coordinates']:
            try:
                polygon = _create_polygon_from_coordinates(data['coordinates'])
                if polygon is None:
                    tls_invalid_count += 1
            except Exception as e:
                print(f"Warning: Could not create polygon for {annotation_id}: {e}")
                tls_invalid_count += 1
        
        all_tls_labels.append(annotation_id)
        all_tls_polygons.append(polygon)
    
    # Process QC annotations
    all_qc_polygons_with_type = []
    qc_invalid_count = 0
    qc_total_count = 0
    
    if qc_annotations:
        for qc_id, qc_data in qc_annotations.items():
            qc_total_count += 1
            qc_poly = None
            if 'coordinates' in qc_data and qc_data['coordinates']:
                try:
                    qc_poly = _create_polygon_from_coordinates(qc_data['coordinates'])
                    if qc_poly is None:
                        qc_invalid_count += 1
                except Exception as e:
                    print(f"Warning: Could not create QC polygon for {qc_id}: {e}")
                    qc_invalid_count += 1
            if qc_poly and not qc_poly.is_empty:
                all_qc_polygons_with_type.append({
                    'polygon': qc_poly,
                    'type_name': qc_data['type']
                })
    
    # Report polygon validation summary
    if tls_invalid_count > 0 or qc_invalid_count > 0:
        print(f"Polygon validation summary:")
        if tls_invalid_count > 0:
            print(f"  - TLS: {tls_invalid_count}/{len(all_tls_polygons)} invalid polygons fixed/skipped")
        if qc_invalid_count > 0:
            print(f"  - QC: {qc_invalid_count}/{qc_total_count} invalid polygons fixed/skipped")
    
    return all_tls_polygons, all_tls_labels, all_qc_polygons_with_type

def calculate_intersection_area(poly1, poly2):
    """
    Calculate the intersection area between two polygons with validation.
    
    Args:
        poly1: First Shapely polygon
        poly2: Second Shapely polygon
        
    Returns:
        float: Intersection area, 0 if no intersection
    """
    try:
        if poly1 is None or poly2 is None or poly1.is_empty or poly2.is_empty:
            return 0.0
        
        # Validate polygons
        if not poly1.is_valid or not poly2.is_valid:
            print("Warning: One or both polygons are invalid for intersection calculation")
            return 0.0
        
        intersection = poly1.intersection(poly2)
        if intersection.is_empty:
            return 0.0
        
        area = intersection.area
        
        # Sanity checks
        if area < 0:
            print(f"Warning: Negative intersection area: {area}")
            return 0.0
        
        # Intersection cannot be larger than either polygon
        if area > poly1.area or area > poly2.area:
            max_possible = min(poly1.area, poly2.area)
            print(f"Warning: Intersection area ({area}) exceeds smaller polygon area ({max_possible})")
            print(f"Polygon 1 area: {poly1.area}, Polygon 2 area: {poly2.area}")
            return max_possible
        
        return area
    except Exception as e:
        print(f"Warning: Error calculating intersection: {e}")
        return 0.0

def calculate_overlap_percentage(polygon_area, intersection_area):
    """
    Calculate what percentage of a polygon is overlapped by another.
    
    Args:
        polygon_area (float): Area of the polygon we're calculating percentage for
        intersection_area (float): Area of the intersection between the polygons
        
    Returns:
        float: Percentage of the polygon that is overlapped (0-100)
    """
    try:
        if polygon_area == 0 or intersection_area == 0:
            return 0.0
        
        # Sanity check: intersection cannot be larger than the polygon
        if intersection_area > polygon_area:
            print(f"Warning: Intersection area ({intersection_area}) > polygon area ({polygon_area})")
            print(f"This suggests invalid geometry. Clamping to 100%.")
            return 100.0
        
        return (intersection_area / polygon_area) * 100
    except Exception as e:
        print(f"Warning: Error calculating overlap percentage: {e}")
        return 0.0

def calculate_tls_qc_intersections(tls_annotations, qc_annotations, sample_name, mpp_x=None, mpp_y=None):
    """
    Calculate intersections between TLS and QC polygons and return detailed statistics.
    
    Args:
        tls_annotations (dict): Dictionary containing TLS annotation data
        qc_annotations (dict): Dictionary containing QC annotation data  
        sample_name (str): Name of the sample being processed
        mpp_x (float, optional): Microns per pixel in X direction
        mpp_y (float, optional): Microns per pixel in Y direction
        
    Returns:
        list: List of dictionaries containing intersection statistics
    """
    print(f"Calculating TLS-QC intersections for sample: {sample_name}")
    if mpp_x is not None and mpp_y is not None:
        print(f"Using MPP values: X={mpp_x:.4f}, Y={mpp_y:.4f} microns/pixel")
    else:
        print("Warning: No MPP values available - areas will be in pixels only")
    start_time = time.time()
    
    # Load and prepare annotations
    all_tls_polygons, all_tls_labels, all_qc_polygons_with_type = _load_and_prepare_annotations(tls_annotations, qc_annotations)
    
    print(f"Loaded {len(all_tls_polygons)} TLS annotations and {len(all_qc_polygons_with_type)} QC annotations")
    
    # Calculate intersections
    intersection_results = []
    
    for tls_idx, tls_polygon in enumerate(all_tls_polygons):
        tls_label = all_tls_labels[tls_idx]
        
        if tls_polygon is None or tls_polygon.is_empty:
            # Record TLS with no valid polygon
            intersection_results.append({
                'sample_name': sample_name,
                'tls_id': tls_label,
                'tls_area_pixels': 0.0,
                'tls_area_microns': 0.0 if mpp_x is not None and mpp_y is not None else None,
                'qc_id': 'N/A',
                'qc_type': 'N/A', 
                'qc_area_pixels': 0.0,
                'qc_area_microns': 0.0 if mpp_x is not None and mpp_y is not None else None,
                'intersection_area_pixels': 0.0,
                'intersection_area_microns': 0.0 if mpp_x is not None and mpp_y is not None else None,
                'tls_overlap_percentage': 0.0,
                'qc_overlap_percentage': 0.0,
                'has_intersection': False
            })
            continue
        
        tls_area_pixels = tls_polygon.area
        tls_area_microns = convert_area_pixels_to_microns(tls_area_pixels, mpp_x, mpp_y)
        has_any_intersection = False
        
        # Check intersections with all QC polygons
        for qc_info in all_qc_polygons_with_type:
            qc_polygon = qc_info['polygon']
            qc_type = qc_info['type_name']
            
            if qc_polygon is None or qc_polygon.is_empty:
                continue
                
            intersection_area_pixels = calculate_intersection_area(tls_polygon, qc_polygon)
            
            if intersection_area_pixels > 0:
                has_any_intersection = True
                qc_area_pixels = qc_polygon.area
                qc_area_microns = convert_area_pixels_to_microns(qc_area_pixels, mpp_x, mpp_y)
                intersection_area_microns = convert_area_pixels_to_microns(intersection_area_pixels, mpp_x, mpp_y)
                
                # Calculate overlap percentages using already-calculated areas
                tls_overlap_pct = calculate_overlap_percentage(tls_area_pixels, intersection_area_pixels)
                qc_overlap_pct = calculate_overlap_percentage(qc_area_pixels, intersection_area_pixels)
                
                # Debug output for problematic cases
                if qc_overlap_pct > 100 or tls_overlap_pct > 100:
                    print(f"DEBUG: Suspicious overlap percentages detected:")
                    print(f"  Sample: {sample_name}")
                    print(f"  TLS ID: {tls_label}")
                    print(f"  QC Type: {qc_type}")
                    print(f"  TLS area: {tls_area_pixels:.2f} pixels²")
                    print(f"  QC area: {qc_area_pixels:.2f} pixels²")
                    print(f"  Intersection area: {intersection_area_pixels:.2f} pixels²")
                    print(f"  TLS overlap: {tls_overlap_pct:.2f}%")
                    print(f"  QC overlap: {qc_overlap_pct:.2f}%")
                
                intersection_results.append({
                    'sample_name': sample_name,
                    'tls_id': tls_label,
                    'tls_area_pixels': tls_area_pixels,
                    'tls_area_microns': tls_area_microns,
                    'qc_id': f"QC_{len(intersection_results)}",  # Generate QC ID
                    'qc_type': qc_type,
                    'qc_area_pixels': qc_area_pixels,
                    'qc_area_microns': qc_area_microns,
                    'intersection_area_pixels': intersection_area_pixels,
                    'intersection_area_microns': intersection_area_microns,
                    'tls_overlap_percentage': tls_overlap_pct,
                    'qc_overlap_percentage': qc_overlap_pct,
                    'has_intersection': True
                })
        
        # If TLS has no intersections with any QC, record it
        if not has_any_intersection:
            intersection_results.append({
                'sample_name': sample_name,
                'tls_id': tls_label,
                'tls_area_pixels': tls_area_pixels,
                'tls_area_microns': tls_area_microns,
                'qc_id': 'None',
                'qc_type': 'None',
                'qc_area_pixels': 0.0,
                'qc_area_microns': 0.0 if mpp_x is not None and mpp_y is not None else None,
                'intersection_area_pixels': 0.0,
                'intersection_area_microns': 0.0 if mpp_x is not None and mpp_y is not None else None,
                'tls_overlap_percentage': 0.0,
                'qc_overlap_percentage': 0.0,
                'has_intersection': False
            })
    
    processing_time = time.time() - start_time
    print(f"Calculated {len(intersection_results)} intersection records in {processing_time:.2f} seconds")
    
    return intersection_results

def save_intersections_to_csv(intersection_results, output_file):
    """
    Save intersection results to a CSV file.
    
    Args:
        intersection_results (list): List of intersection dictionaries
        output_file (str): Path to output CSV file
    """
    print(f"Saving intersection results to: {output_file}")
    
    # Define CSV headers with both pixel and micron measurements
    headers = [
        'sample_name',
        'tls_id', 
        'tls_area_pixels',
        'tls_area_microns',
        'qc_id',
        'qc_type',
        'qc_area_pixels',
        'qc_area_microns',
        'intersection_area_pixels',
        'intersection_area_microns',
        'tls_overlap_percentage',
        'qc_overlap_percentage', 
        'has_intersection'
    ]
    
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        writer.writerows(intersection_results)
    
    print(f"Successfully saved {len(intersection_results)} records to {output_file}")



if __name__ == "__main__":
    # Check if we have the right number of arguments
    if len(sys.argv) < 4:
        print("Usage: python calc_intersect.py <tls_json_file> <qc_json_file> <output_csv_file> [wsi_file]")
        print("Example: python calc_intersect.py tls.geojson qc.geojson intersections.csv sample.svs")
        print("Note: WSI file is optional but recommended for MPP normalization")
        sys.exit(1)
    
    # Get arguments
    tls_json_file = sys.argv[1]
    qc_json_file = sys.argv[2]
    output_csv_file = sys.argv[3]
    wsi_file = sys.argv[4] if len(sys.argv) > 4 else None
    
    # Extract sample name from file path
    sample_name = os.path.splitext(os.path.basename(tls_json_file))[0]
    print(f"Processing sample: {sample_name}")
    
    # Make sure output directory exists
    os.makedirs(os.path.dirname(output_csv_file), exist_ok=True)
    
    try:
        # Extract MPP values from WSI if provided
        mpp_x, mpp_y = None, None
        if wsi_file:
            if os.path.exists(wsi_file):
                mpp_x, mpp_y = extract_mpp_from_wsi(wsi_file)
            else:
                print(f"Warning: WSI file not found: {wsi_file}")
        else:
            print("No WSI file provided - areas will be reported in pixels only")
        
        # Process TLS file
        print(f"Loading TLS file: {tls_json_file}")
        with open(tls_json_file, 'r') as file:
            tls_json = json.load(file)
        tls_annotations = tls_annotations_json(tls_json)
        print(f"Loaded TLS annotations: {len(tls_annotations)} items")
        
        # Process QC file
        print(f"Loading QC file: {qc_json_file}")
        with open(qc_json_file, 'r') as file:
            qc_json = json.load(file)
        qc_annotations = qc_annotations_json(qc_json)
        print(f"Loaded QC annotations: {len(qc_annotations)} items")
        
        # Calculate intersections with MPP normalization
        intersection_results = calculate_tls_qc_intersections(
            tls_annotations, qc_annotations, sample_name, mpp_x, mpp_y
        )
        
        # Save to CSV
        save_intersections_to_csv(intersection_results, output_csv_file)
        
        print(f"\nProcessing complete!")
        print(f"Sample: {sample_name}")
        print(f"TLS annotations: {len(tls_annotations)}")
        print(f"QC annotations: {len(qc_annotations)}")
        print(f"Intersection records: {len(intersection_results)}")
        if mpp_x is not None and mpp_y is not None:
            print(f"MPP normalization: X={mpp_x:.4f}, Y={mpp_y:.4f} microns/pixel")
        else:
            print("MPP normalization: Not available")
        print(f"Output saved to: {output_csv_file}")
        
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON format - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)