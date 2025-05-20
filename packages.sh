# Created a new conda environment called grandqc
conda create -n grandqc python==3.10 -y
conda activate grandqc

# Downloading GrandQC and installing requirements
git clone https://github.com/cpath-ukk/grandqc.git && cd grandqc
pip install -r requirements.txt
conda install openslide openslide-python

#making the Snakemake workflow
conda install snakemake