#!/bin/bash
#SBATCH --job-name=grandqc_run_tis
#SBATCH --output=grandqc_run_tis.out
#SBATCH --error=grandqc_run_tis.err
#SBATCH --time=48:00:00
#SBATCH --qos=normal
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1


# Load conda and activate the environment
source /ddn_exa/campbell/pgupta/miniforge3/etc/profile.d/conda.sh
conda activate grandqc

echo Successfully activated grandqc conda environment

# Run the tissue detection script
cd /ddn_exa/campbell/pgupta/grandqc/01_WSI_inference_OPENSLIDE_QC/
# python wsi_tis_detect.py \
#     --slide_folder /ddn_exa/campbell/share/datasets/temp_TCGA/ \
#     --output_dir /ddn_exa/campbell/pgupta/2025-HistoPath-TLS-Purav/TCGA_out_2/ > /ddn_exa/campbell/pgupta/2025-HistoPath-TLS-Purav/grandqc_log.txt


python main.py \
    --slide_folder /ddn_exa/campbell/share/datasets/TCGA-histopath/ \
    --output_dir /ddn_exa/campbell/pgupta/2025-HistoPath-TLS-Purav/TCGA_out/ \
    --create_geojson "Y" \
    --mpp_model "1.5" > /ddn_exa/campbell/pgupta/2025-HistoPath-TLS-Purav/TCGA_out/grandqc_log.txt