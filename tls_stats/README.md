# TLS-QC Intersection Analysis with MPP Normalization

This directory contains tools for calculating intersections between Tertiary Lymphoid Structures (TLS) and Quality Control (QC) annotations, with support for MPP (Microns Per Pixel) normalization.

## Features

- **Intersection Analysis**: Calculate overlaps between TLS and QC polygons
- **MPP Normalization**: Convert pixel-based measurements to physical units (microns)
- **Dual Output**: Results include both pixel and micron measurements
- **Batch Processing**: Process multiple samples automatically
- **Quality Control**: Support for various QC annotation types (Fold, Darkspot, PenMarking, etc.)

## Quick Start

### 1. Configure WSI File Path
Edit `run_intersection_analysis.sh` and update the `WSI_DIR` variable to point to your whole slide image files:

```bash
WSI_DIR="/path/to/your/wsi_files"  # Update this path
```

### 2. Run Batch Analysis
```bash
cd /ddn_exa/campbell/pgupta/2025-HistoPath-TLS-Purav/tls_stats
./run_intersection_analysis.sh
```

### 3. Run Single Sample
```bash
python calc_intersect.py tls.geojson qc.geojson output.csv [wsi_file.svs]
```

## File Structure

### Expected Input Files
- **TLS files**: `mod_{sample_name}.geojson` (in `../mod_TLS_json/`)
- **QC files**: `{sample_name}.svs.geojson` (in `../qc_geojson/`)
- **WSI files**: `{sample_name}.svs` (in your configured WSI directory)

### Output Files
- **Individual results**: `{sample_name}_intersections.csv`
- **Combined results**: `all_samples_intersections.csv`

## Output Format

The CSV output includes the following columns:

### Identification
- `sample_name`: Name of the processed sample
- `tls_id`: TLS annotation identifier
- `qc_id`: QC annotation identifier
- `qc_type`: Type of QC issue (Fold, Darkspot, etc.)

### Area Measurements (Dual Units)
- `tls_area_pixels` / `tls_area_microns`: TLS annotation area
- `qc_area_pixels` / `qc_area_microns`: QC annotation area  
- `intersection_area_pixels` / `intersection_area_microns`: Overlap area

### Overlap Statistics
- `tls_overlap_percentage`: % of TLS covered by QC
- `qc_overlap_percentage`: % of QC covered by TLS
- `has_intersection`: Boolean indicating if overlap exists

## MPP Normalization

### How It Works
The script extracts MPP values from WSI metadata using OpenSlide:
- Tries multiple metadata formats (Aperio, OpenSlide, TIFF)
- Converts pixel areas to micron areas using: `area_microns = area_pixels × mpp_x × mpp_y`

### Graceful Degradation
- **With WSI files**: Areas reported in both pixels and microns
- **Without WSI files**: Areas reported in pixels only (with warnings)
- **Invalid MPP**: Falls back to pixel-only measurements

### WSI Format Support
Supports major WSI formats through OpenSlide:
- Aperio (.svs)
- Hamamatsu (.ndpi)
- Leica (.scn)
- And others supported by OpenSlide

## Error Handling

### Common Issues
1. **WSI directory not found**: Updates required in shell script
2. **Missing WSI files**: Script continues with pixel-only measurements
3. **Invalid MPP metadata**: Falls back gracefully with warnings
4. **Missing annotation files**: Skips problematic samples with warnings
5. **Invalid polygon geometry**: Common with annotation data - automatically fixed when possible

### Polygon Validation Warnings
The script automatically detects and fixes common polygon issues:
- **Self-intersecting polygons**: Fixed using Shapely's buffer(0) operation
- **Zero-area polygons**: Skipped (often from degenerate annotations)
- **Invalid coordinates**: Gracefully handled and skipped

These warnings are **normal** and indicate the script is working correctly. The validation summary shows how many polygons were processed vs. fixed/skipped.

### Troubleshooting
- Check file paths and permissions
- Verify WSI files are readable by OpenSlide
- Ensure annotation files are valid GeoJSON format
- Review console output for specific error messages

## Dependencies

- Python packages: `openslide-python`, `shapely`, `numpy`
- System: OpenSlide library
- Input: GeoJSON annotation files, WSI files

## Example Usage

```bash
# Process single sample with MPP normalization
python calc_intersect.py \
    ../mod_TLS_json/mod_sample1.geojson \
    ../qc_geojson/sample1.svs.geojson \
    sample1_results.csv \
    /path/to/wsi/sample1.svs

# Process all samples
./run_intersection_analysis.sh
```

## Output Example

```csv
sample_name,tls_id,tls_area_pixels,tls_area_microns,qc_id,qc_type,intersection_area_pixels,intersection_area_microns,tls_overlap_percentage,has_intersection
sample1,TLS_001,178134904.0,43821.2,QC_001,Fold,50234.5,12.3,0.028,true
```

For questions or issues, check the console output for detailed error messages and warnings.
