import os
import glob
from os.path import basename, splitext, join

# Load config
configfile: "config.yaml"

RESULTS_DIR = config["results_dir"]          # Pull from config.yaml for portability
WSI_DIR = config['input_dir']                # Directory with WSI images
TLS_GEOJSON_DIR = config['tls_geojson_dir']  # Directory with TLS geojson files
QC_GEOJSON_DIR = config['qc_geojson_dir']    # Directory with QC geojson files
WSI_LEVEL = config.get("wsi_level", 0)       # Default to 0 if not specified in config

# Dynamically find all WSI files
WSI_FILES = glob.glob(os.path.join(WSI_DIR, "*.svs"))
print(f"Found {len(WSI_FILES)} WSI files: {[basename(f) for f in WSI_FILES]}")

# Extract base names (without .svs extension) for all WSI files
SAMPLE_IDS = [splitext(basename(f))[0] for f in WSI_FILES]
print(f"Sample IDs: {SAMPLE_IDS}")

# Wildcard constraints
wildcard_constraints:
    sample = r"[^/]+",  # Sample names cannot contain forward slashes
    level = r"\d+"      # Level must be digits only

# Create necessary directories
os.makedirs(join(RESULTS_DIR, "logs"), exist_ok=True)
os.makedirs(join(RESULTS_DIR, "tls_images"), exist_ok=True)
os.makedirs(join(RESULTS_DIR, "tls_qc_images"), exist_ok=True)


rule all:
    input:
        # Generate all TLS-only images for all samples
        expand(os.path.join(RESULTS_DIR, "tls_images", "{sample}_level{level}_TLS.png"), 
               sample=SAMPLE_IDS, level=WSI_LEVEL),
        # Generate all TLS+QC images for all samples
        expand(os.path.join(RESULTS_DIR, "tls_qc_images", "{sample}_level{level}_TLS_QC.png"), 
               sample=SAMPLE_IDS, level=WSI_LEVEL),
        # Completion flag
        os.path.join(RESULTS_DIR, "visualization.done")

rule generate_tls_image:
    input:
        tls_geojson = os.path.join(TLS_GEOJSON_DIR, "{sample}.geojson"),
        wsi = os.path.join(WSI_DIR, "{sample}.svs")
    output:
        image = os.path.join(RESULTS_DIR, "tls_images", "{sample}_level{level}_TLS.png")
    resources:
        mem_mb=16000,
        runtime=180,
        slurm_partition="guest"
    log:
        os.path.join(RESULTS_DIR, "logs", "{sample}_level{level}_tls.log")
    shell:
        """
        python main_visualizer_tls.py {input.tls_geojson} {input.wsi} {output.image} {wildcards.level} > {log} 2>&1
        """

rule generate_tls_qc_image:
    input:
        tls_geojson = os.path.join(TLS_GEOJSON_DIR, "{sample}.geojson"),
        qc_geojson = os.path.join(QC_GEOJSON_DIR, "{sample}.svs.geojson"),
        wsi = os.path.join(WSI_DIR, "{sample}.svs")
    output:
        image = os.path.join(RESULTS_DIR, "tls_qc_images", "{sample}_level{level}_TLS_QC.png")
    resources:
        mem_mb=32000,
        runtime=180,
        slurm_partition="guest"
    log:
        os.path.join(RESULTS_DIR, "logs", "{sample}_level{level}_tls_qc.log")
    shell:
        """
        python main_visualizer_tls_qc.py {input.tls_geojson} {input.qc_geojson} {input.wsi} {output.image} {wildcards.level} > {log} 2>&1
        """

rule finalize_visualization:
    input:
        tls_images = expand(os.path.join(RESULTS_DIR, "tls_images", "{sample}_level{level}_TLS.png"), 
                           sample=SAMPLE_IDS, level=WSI_LEVEL),
        tls_qc_images = expand(os.path.join(RESULTS_DIR, "tls_qc_images", "{sample}_level{level}_TLS_QC.png"), 
                              sample=SAMPLE_IDS, level=WSI_LEVEL)
    output:
        touch(os.path.join(RESULTS_DIR, "visualization.done"))
    resources:
        mem_mb=4000,
        runtime=30,
        slurm_partition="guest"
    log:
        os.path.join(RESULTS_DIR, "logs", "finalization.log")
    shell:
        """
        echo "Visualization completed successfully at $(date)" > {log}
        echo "Generated $(ls {RESULTS_DIR}/tls_images/*.png | wc -l) TLS images" >> {log}
        echo "Generated $(ls {RESULTS_DIR}/tls_qc_images/*.png | wc -l) TLS+QC images" >> {log}
        echo "Total images: $(ls {RESULTS_DIR}/tls_images/*.png {RESULTS_DIR}/tls_qc_images/*.png | wc -l)" >> {log}
        echo "WSI level used: {WSI_LEVEL}" >> {log}
        echo "Processed samples: {SAMPLE_IDS}" >> {log}
        """