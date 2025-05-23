#!/bin/bash
# Script to run the TLS visualization Snakemake workflow
rm TLS_gallery/visualization.done

# Load any necessary modules or environment settings
echo "Starting TLS visualization workflow"

# Run Snakemake with the profile
snakemake --profile profiles/default

# Check if the workflow completed successfully
if [ $? -eq 0 ]; then
    echo "Visualization workflow completed successfully!"
    echo "Results are available in: $(grep 'results_dir' config.yaml | cut -d':' -f2 | tr -d ' "')"
else
    echo "Visualization workflow encountered an error. Please check the logs."
    echo "Log files are in: $(grep 'results_dir' config.yaml | cut -d':' -f2 | tr -d ' "')/logs/"
fi
