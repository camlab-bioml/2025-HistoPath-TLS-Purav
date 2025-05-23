import json
import time
import sys
import os
from visualize_annotations_tls_qc import display_annotations_tls_qc
from parse_json import tls_annotations_json, qc_annotations_json, modify_json

def generate_tls_qc_visualization(tls_json_file, qc_json_file, image_path, output_file, wsi_level=0):
    """
    Generate TLS+QC visualization and save it to the specified output file.
    
    Args:
        tls_json_file (str): Path to the TLS geojson file
        qc_json_file (str): Path to the QC geojson file
        image_path (str): Path to the WSI image file
        output_file (str): Path to save the output visualization
        wsi_level (int): WSI level for visualization (0 is highest resolution)
    """
    print(f"Starting TLS+QC visualization...")
    start_time = time.time()
    
    # Process TLS file
    print(f"Processing TLS file: {tls_json_file}")
    modified_json = modify_json(tls_json_file)

    # Process QC file
    print(f"Processing QC file: {qc_json_file}")
    with open(qc_json_file, 'r') as file:
        qc_json = json.load(file)
    
    qc_annotations = qc_annotations_json(qc_json)
    print(f"Loaded QC annotations: {len(qc_annotations)} items")

    # Load TLS annotations
    annotations = tls_annotations_json(modified_json)
    print(f"Loaded TLS annotations: {len(annotations)} items")
    
    # Generate TLS+QC visualization
    print(f"Generating TLS+QC visualization at level {wsi_level}...")
    
    display_annotations_tls_qc(
        qc_annotations=qc_annotations,
        tls_annotations=annotations,
        image_path=image_path,
        save_path=output_file,
        with_mask=True,
        num_per_row=7,
        file_format="png",
        level=wsi_level
    )
    
    print(f"TLS+QC visualization complete. Image saved to: {output_file}")
    print(f"Processing time: {time.time() - start_time:.2f} seconds")

if __name__ == "__main__":
    # Check if we have the right number of arguments
    if len(sys.argv) < 5:
        print("Usage: python main_visualizer_tls_qc.py <tls_json_file> <qc_json_file> <wsi_file> <output_file> [wsi_level]")
        sys.exit(1)
    
    # Get arguments
    tls_json_file = sys.argv[1]
    qc_json_file = sys.argv[2]
    wsi_file = sys.argv[3]
    output_file = sys.argv[4]
    
    # Optional wsi_level argument, default to 0
    wsi_level = 0
    if len(sys.argv) > 5:
        try:
            wsi_level = int(sys.argv[5])
        except ValueError:
            print(f"Warning: Invalid WSI level '{sys.argv[5]}', using default level 0")
    
    # Make sure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Generate visualization
    generate_tls_qc_visualization(tls_json_file, qc_json_file, wsi_file, output_file, wsi_level)
