#!/bin/bash

# Script to run TLS visualization pipeline on all WSI samples
# This script provides different strategies for running Snakemake jobs

# Set working directory
cd /ddn_exa/campbell/pgupta/2025-HistoPath-TLS-Purav

echo "=== TLS Visualization Pipeline - All Samples ==="
echo "Working directory: $(pwd)"
echo "Date: $(date)"

# Check if WSI directory exists and count files
WSI_DIR=$(grep "input_dir:" config.yaml | cut -d'"' -f2)
echo "WSI directory: $WSI_DIR"

if [ -d "$WSI_DIR" ]; then
    WSI_COUNT=$(ls -1 "$WSI_DIR"/*.svs 2>/dev/null | wc -l)
    echo "Found $WSI_COUNT WSI files in $WSI_DIR"
    if [ $WSI_COUNT -eq 0 ]; then
        echo "ERROR: No .svs files found in WSI directory!"
        exit 1
    fi
else
    echo "ERROR: WSI directory does not exist: $WSI_DIR"
    exit 1
fi

# Function to show usage
show_usage() {
    echo ""
    echo "Usage: $0 [strategy]"
    echo ""
    echo "Strategies:"
    echo "  1. dry-run     - Show what jobs would be submitted (recommended first)"
    echo "  2. sequential  - Run one job at a time (safe, slow)"
    echo "  3. parallel    - Run multiple jobs in parallel (fast, needs resources)"
    echo "  4. custom      - Custom job limits and resources"
    echo ""
    echo "Examples:"
    echo "  $0 dry-run     # See what will be executed"
    echo "  $0 sequential  # Run jobs one by one"
    echo "  $0 parallel    # Run up to 10 jobs in parallel"
    echo ""
}

# Parse command line arguments
STRATEGY=${1:-"help"}

case $STRATEGY in
    "dry-run"|"dryrun"|"dry")
        echo ""
        echo "=== DRY RUN - Showing what would be executed ==="
        snakemake --profile profiles/default --dry-run --printshellcmds
        echo ""
        echo "To actually run the pipeline, use one of the execution strategies."
        show_usage
        ;;
        
    "sequential"|"seq")
        echo ""
        echo "=== SEQUENTIAL EXECUTION ==="
        echo "Running jobs one at a time..."
        snakemake --profile profiles/default --jobs 1 --printshellcmds
        ;;
        
    "parallel"|"par")
        echo ""
        echo "=== PARALLEL EXECUTION ==="
        echo "Running up to 100 jobs in parallel (from profile config)..."
        snakemake --profile profiles/default --printshellcmds
        ;;
        
    "custom")
        echo ""
        echo "=== CUSTOM EXECUTION ==="
        echo "Enter custom parameters:"
        read -p "Number of parallel jobs (default 5): " JOBS
        JOBS=${JOBS:-5}
        
        read -p "Additional snakemake options (optional): " EXTRA_OPTS
        
        echo "Running with $JOBS parallel jobs..."
        snakemake --profile profiles/default --jobs $JOBS --printshellcmds $EXTRA_OPTS
        ;;
        
    "help"|"-h"|"--help"|*)
        echo ""
        echo "This script will process ALL WSI files in the input directory."
        echo "Each WSI file will be processed as a separate Slurm job."
        show_usage
        exit 0
        ;;
esac

# Check if execution was successful
if [ $? -eq 0 ]; then
    echo ""
    echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
    echo "Check results in: $(grep "results_dir:" config.yaml | cut -d'"' -f2)"
    echo "TLS images: $(grep "results_dir:" config.yaml | cut -d'"' -f2)/tls_images/"
    echo "TLS+QC images: $(grep "results_dir:" config.yaml | cut -d'"' -f2)/tls_qc_images/"
else
    echo ""
    echo "=== PIPELINE FAILED ==="
    echo "Check the logs in: $(grep "results_dir:" config.yaml | cut -d'"' -f2)/logs/"
    exit 1
fi
