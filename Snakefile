import os

# Load config
configfile: "config.yaml"

SVS_DIR = config["input_dir"]
RESULTS_DIR = config["results_dir"]
GRANDQC_PATH = config["grandqc_path"]

rule all:
    input:
        os.path.join(RESULTS_DIR, "artifact_segmentation.done")

rule run_tis:
    output:
        touch(os.path.join(RESULTS_DIR, "tissue_segmentation.done"))
    shell:
        """
        cd {GRANDQC_PATH} && \
        python wsi_tis_detect.py --slide_folder {SVS_DIR} --output_dir {RESULTS_DIR}
        """

rule run_art:
    input:
        os.path.join(RESULTS_DIR, "tissue_segmentation.done")
    output:
        touch(os.path.join(RESULTS_DIR, "artifact_segmentation.done"))
    shell:
        """
        cd {GRANDQC_PATH} && \
        python main.py --slide_folder {SVS_DIR} --output_dir {RESULTS_DIR} --create_geojson "Y" --mpp_model "1.5"
        """