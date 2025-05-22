import os

# Load config
configfile: "config.yaml"

SVS_DIR = config["input_dir"]
RESULTS_DIR = config["results_dir"]
GRANDQC_PATH = config["grandqc_path"]

# Create log directory
LOG_DIR = os.path.join(RESULTS_DIR, "logs")
os.makedirs(LOG_DIR, exist_ok=True)

# Make directory for tissue masks if it doesn't exist
TIS_DET_MASK_DIR = os.path.join(RESULTS_DIR, "tis_det_mask")
os.makedirs(TIS_DET_MASK_DIR, exist_ok=True)

rule all:
    input:
        os.path.join(RESULTS_DIR, "artifact_segmentation.done")

rule run_tis:
    output:
        touch(os.path.join(RESULTS_DIR, "tissue_segmentation.done"))
    resources:
        mem_mb=16000,
        runtime=360,
        gpus=1,
        slurm_partition="gpu",
        nodelist="gpu1" 
    log:
        os.path.join(LOG_DIR, "tissue_segmentation.log")
    shell:
        """
        # Activate conda environment manually to ensure all paths are correct
        source /ddn_exa/campbell/pgupta/miniforge3/etc/profile.d/conda.sh && 
        conda activate grandqc &&
        cd {GRANDQC_PATH} && 
        python wsi_tis_detect.py --slide_folder {SVS_DIR} --output_dir {RESULTS_DIR} > {log} 2>&1
        """

rule run_art:
    input:
        os.path.join(RESULTS_DIR, "tissue_segmentation.done")
    output:
        touch(os.path.join(RESULTS_DIR, "artifact_segmentation.done"))
    resources:
        mem_mb=16000,
        runtime=720,
        gpus=1,
        slurm_partition="gpu",
        nodelist="gpu1"
    log:
        os.path.join(LOG_DIR, "artifact_segmentation.log")
    shell:
        """
        # Activate conda environment manually to ensure all paths are correct
        source /ddn_exa/campbell/pgupta/miniforge3/etc/profile.d/conda.sh && 
        conda activate grandqc &&
        cd {GRANDQC_PATH} && 
        python main.py --slide_folder {SVS_DIR} --output_dir {RESULTS_DIR} --create_geojson "Y" --mpp_model "1.5" > {log} 2>&1
        """