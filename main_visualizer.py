import json
import cProfile
from visualize_annotations import tls_annotations_json, qc_annotations_json, display_annotations
from mod_json import modify_json

# Example usage
if __name__ == "__main__":
    # Replace with your JSON file path
    tls_json_file = "TCGA-IB-AAUU-01Z-00-DX1.83137319-8229-4682-8BD3-9B9A5C6C997A.hooknet.geojson"
    
    # Image path
    image_path = '/Users/puravgupta/Documents/Research/Pancreas_TLS_ML/WSI/TCGA-IB-AAUU-01Z-00-DX1.83137319-8229-4682-8BD3-9B9A5C6C997A.svs' 
    
    modified_json = modify_json(tls_json_file)

    qc_json_file = "TCGA-IB-AAUU-01Z-00-DX1.83137319-8229-4682-8BD3-9B9A5C6C997A.svs.geojson"

    with open(qc_json_file, 'r') as file:
        qc_json = json.load(file)
    
    qc_annotations = qc_annotations_json(qc_json)

    # Load annotations
    annotations = tls_annotations_json(modified_json)
    
    # Set the WSI level (0 is highest resolution, higher numbers mean lower resolution)
    # Try different levels to see performance impact
    wsi_level = 0
    
    # Use cProfile to analyze performance
    profiler = cProfile.Profile()
    profiler.enable()
    
    # Display annotations
    display_annotations(
        qc_annotations=qc_annotations,
        tls_annotations=annotations,
        image_path=image_path,
        save_path=f"{image_path.split('/')[-1].split('.')[0]}_level{wsi_level}.png",
        with_mask=True,
        num_per_row=7,
        file_format="png",
        level=wsi_level
    )
    
    profiler.disable()
    profiler.print_stats(sort='cumtime')
