import json
import time
import sys
import os
from visualize_annotations_tls import display_annotations_tls
from parse_json import tls_annotations_json, modify_json

def generate_tls_visualization(tls_json_file, image_path, output_file, wsi_level=0):
    """
    Generate TLS-only visualization and save it to the specified output file.
    
    Args:
        tls_json_file (str): Path to the TLS geojson file
        image_path (str): Path to the WSI image file
        output_file (str): Path to save the output visualization
        wsi_level (int): WSI level for visualization (0 is highest resolution)
    """
    print(f"Starting TLS-only visualization...")
    start_time = time.time()
    
    # Process TLS file
    print(f"Processing TLS file: {tls_json_file}")
    modified_json = modify_json(tls_json_file)

    # Load annotations
    annotations = tls_annotations_json(modified_json)
    print(f"Loaded TLS annotations: {len(annotations)} items")
    
    # Generate TLS-only visualization
    print(f"Generating TLS visualization at level {wsi_level}...")
    
    display_annotations_tls(
        tls_annotations=annotations,
        image_path=image_path,
        save_path=output_file,
        with_mask=True,
        num_per_row=7,
        file_format="png",
        level=wsi_level
    )
    
    print(f"TLS visualization complete. Image saved to: {output_file}")
    print(f"Processing time: {time.time() - start_time:.2f} seconds")

if __name__ == "__main__":
    # Check if we have the right number of arguments
    if len(sys.argv) < 4:
        print("Usage: python main_visualizer_tls.py <tls_json_file> <wsi_file> <output_file> [wsi_level]")
        sys.exit(1)
    
    # Get arguments
    tls_json_file = sys.argv[1]
    wsi_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Optional wsi_level argument, default to 0
    wsi_level = 0
    if len(sys.argv) > 4:
        try:
            wsi_level = int(sys.argv[4])
        except ValueError:
            print(f"Warning: Invalid WSI level '{sys.argv[4]}', using default level 0")
    
    # Make sure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Generate visualization
    generate_tls_visualization(tls_json_file, wsi_file, output_file, wsi_level)
