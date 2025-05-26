#!/bin/bash

echo "Starting Snakemake workflow..."
snakemake --profile profiles/default

echo "Snakemake workflow completed."

OUTPUT_DIR="/ddn_exa/campbell/pgupta/2025-HistoPath-TLS-Purav/tls_stats/intersection_results"

# Combine all results into a single file
COMBINED_CSV="$OUTPUT_DIR/all_samples_intersections.csv"
echo ""
echo "Combining all results into: $COMBINED_CSV"

# Create combined CSV with header
first_file=true
for csv_file in "$OUTPUT_DIR"/*.csv; do
    if [ ! -f "$csv_file" ]; then
        echo "No CSV files found in $OUTPUT_DIR"
        continue
    fi
    
    # Skip the combined file itself if it exists
    if [[ "$csv_file" == "$COMBINED_CSV" ]]; then
        continue
    fi
    
    if [ "$first_file" = true ]; then
        # Include header from first file
        cat "$csv_file" > "$COMBINED_CSV"
        first_file=false
    else
        # Skip header for subsequent files
        tail -n +2 "$csv_file" >> "$COMBINED_CSV"
    fi
done

if [ -f "$COMBINED_CSV" ]; then
    TOTAL_RECORDS=$(tail -n +2 "$COMBINED_CSV" | wc -l)
    echo "Combined CSV created with $TOTAL_RECORDS total records"
    echo "Output saved to: $COMBINED_CSV"
else
    echo "No CSV files were found to combine"
fi

echo "Analysis complete!"