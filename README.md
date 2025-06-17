# README: Moment Maps Generation from ALMA Data using CASA

## Overview

This repository provides a script for generating moment maps (Moment 0, 1, and 2) for various molecular lines (CO J=0-1, 2-1, 3-2, 4-3, HCO+, CI, C+) from ALMA data. The scripts utilize the Common Astronomy Software Applications (CASA) package to process the data and create maps for the Circum-Nuclear Disk (CND) scale.

## Prerequisites

- **CASA**: Ensure you have CASA installed on your system. Download it from the [CASA website].
- **Python**: Python 3.x is recommended for running the scripts.
- **ALMA Data**: Obtain the necessary ALMA datasets in FITS format.

## Usage

### Script Configuration

1. **Edit the Script**: The main script is configured to import FITS files, process them, and generate moment maps. You can adjust the following parameters:
   - `snaps`: List of snapshot identifiers to process.
   - `incl`: Inclination angle used in the processing.
   - `pixel`: Size of the image in pixels.

### Running the Script

To run the script, execute the following command in the CASA environment:

```bash
casa --nologger --nogui -c CASA_script.py
```

### Workflow Steps

1. **Import FITS Files**: The script imports FITS files into CASA image format.
2. **Get Beam Size**: Retrieves beam size from the header of the images.
3. **Set Header Information**: Updates header information to define the units (Jy/pixel).
4. **Generate Moment Maps**: Creates moment maps and exports them in FITS format.
5. **Clean Up**: Removes temporary files generated during processing.

### Example Script Structure

Here’s an overview of key sections in the script:

```python
import numpy as np
import astropy.units as u

# Configuration
pixel = 1024
snaps = ['030', '040', '050', '060', '070', '080']
incl = 41

# Import FITS files into CASA
for snap in snaps:
    FIT_IM = [
        f'image_s{snap}_line1_CI_i{incl}_pix{pixel}_Jypixel.fits',
        f'image_s{snap}_line1_CO_i{incl}_pix{pixel}_Jypixel.fits',
    ]
    for fitsim in FIT_IM:
        importfits(fitsimage=fitsim, imagename=fitsim.split('.')[0] + '.image')

# Generate Moment Maps
for image in IMG_IM:
    immoments(imagename=image, moments=[0, 1, 2], outfile=image.split('.')[0])

# Export to FITS format
for ext in ext_casa:
    exportfits(imagename=image, fitsimage=image.split('.')[0] + f'_{ext}.fits')
```

## Output

The script generates the following moment maps:

- **Moment 0**: Integrated intensity map.
- **Moment 1**: Velocity field map.
- **Moment 2**: Velocity dispersion map.

These maps are exported in FITS format for further analysis.

## Cleaning Up

The script includes commands to remove temporary files to keep the workspace tidy:

```python
os.system('rm -rf *.last')
os.system('rm -rf *.image *.integrated *.weighted_coord *.weighted_dispersion_coord')
```


**Radio Astronomy Data Processing Script**

The script automates the processing of radio astronomy data to produce scientifically usable outputs, such as cleaned images or statistical summaries. It demonstrates a complete workflow, including:

Importing and inspecting raw data
Flagging corrupted or noisy data
Calibrating instrumental and environmental effects
Creating sky images using the CLEAN algorithm
Extracting basic image statistics
Organizing outputs for further analysis

Requirements

CASA: Version 6.x or later, available from the National Radio Astronomy Observatory (NRAO).
Python: Compatible with CASA’s Python environment (typically Python 2.7 or 3.x, depending on CASA version).
Operating System: Linux, macOS, or Windows (with CASA installed).
Dependencies: CASA modules (casatasks, casatools), Python standard libraries (os, glob).
Disk Space: Sufficient space for input Measurement Set (MS) and output files (e.g., calibration tables, FITS images).

Installation

Install CASA following the instructions on the NRAO website.
Ensure the CASA Python environment is accessible (e.g., run casapy or activate CASA’s Python).
Place the script (casa_data_processing.py) in your working directory.
Create a data/ directory for the input Measurement Set and an outputs/ directory for results (the script creates outputs/ if it doesn’t exist).

Usage

Prepare Input Data:

For demonstration, the script assumes a dummy Measurement Set (data/dummy_data.ms).
For real data, replace with an actual Measurement Set and uncomment the importasdm task to convert raw telescope data (e.g., ALMA ASDM format) to MS format.


Run the Script:
casapy -c casa_data_processing.py

Alternatively, run in a Python environment with CASA modules:
python casa_data_processing.py


Check Outputs:Outputs are saved in the outputs/ directory. See the "Outputs" section below for details.


File Structure
project_directory/
├── data/
│   └── dummy_data.ms          # Input Measurement Set (dummy or real data)
├── outputs/                   # Output directory (created by script)
│   ├── obs_summary.txt        # Observation summary
│   ├── amp_vs_time.png        # Amplitude vs. time plot
│   ├── bandpass.cal           # Bandpass calibration table
│   ├── gain.cal               # Gain calibration table
│   ├── target_image.*         # Cleaned image files
│   ├── target_image.fits      # Exported FITS image
│   └── image_stats.txt        # Image statistics
└── casa_data_processing.py    # Main script

Customization
To adapt the script for real data, modify the following parameters based on your dataset:

File Paths:

Update data_dir and ms_file to point to your Measurement Set.
Example: ms_file = './data/my_observation.ms'


Observation Parameters:

Use listobs output to identify field names (e.g., CAL_BP, CAL_PHASE, TARGET) and spectral windows (spw).
Update field and spw in tasks like bandpass, gaincal, applycal, and tclean.


Imaging Parameters:

Adjust imsize (image size in pixels), cell (pixel size in arcseconds), and threshold (cleaning threshold) in tclean based on your data’s resolution and sensitivity.
Example: For high-resolution ALMA data, set cell='0.1arcsec' and imsize=[512, 512].


Calibration:

Replace CAL_BP, CAL_PHASE, and TARGET with actual source names from listobs.
Adjust refant (reference antenna) and solint (solution interval) based on data quality.


Spectral Analysis (Optional):

Uncomment the spectral cube section in the script to create and analyze a spectral cube (e.g., for HI or CO lines).
Update region in imval for specific spectral extractions.



Outputs
The script generates the following files in the outputs/ directory:

obs_summary.txt: Text file summarizing observation details (antennas, sources, spectral windows).
amp_vs_time.png: Plot of amplitude vs. time for data quality inspection.
bandpass.cal: Calibration table for frequency-dependent corrections.
gain.cal: Calibration table for time-dependent amplitude and phase corrections.
target_image.*: CASA image files (e.g., .image, .residual) from tclean.
target_image.fits: Exported FITS image for external analysis.
image_stats.txt: Statistics of the cleaned image (peak flux, RMS noise).

Example Workflow

Import and Inspect:

Imports a Measurement Set and generates a summary (obs_summary.txt).
Creates a diagnostic plot (amp_vs_time.png) to check for RFI or anomalies.


Flagging:

Flags shadowed data, high-amplitude outliers, and a specific bad antenna (e.g., ANT1).


Calibration:

Applies bandpass calibration using a bright calibrator (CAL_BP).
Performs gain calibration for amplitude and phase (CAL_PHASE).
Applies calibrations to the science target (TARGET).


Imaging:

Creates a cleaned image using the CLEAN algorithm (target_image).
Exports the image as a FITS file (target_image.fits).


Analysis:

Computes image statistics (peak flux, RMS noise) and saves them (image_stats.txt).


Cleanup:

Lists all generated files for reference.
Optional cleanup of temporary calibration tables (commented out).



Troubleshooting

Missing MS File: If dummy_data.ms is not found, replace with a real Measurement Set or use importasdm for raw data.
Calibration Failures: Check listobs output for correct calibrator names and ensure sufficient signal-to-noise ratio (minsnr).
Imaging Artifacts: Adjust tclean parameters (e.g., niter, threshold, robust) or re-flag data if sidelobes persist.
CASA Errors: Ensure CASA is installed correctly and run in a compatible environment. Check CASA logs for detailed error messages.

Notes

The script is designed for demonstration with a dummy dataset. For real data, consult telescope documentation (e.g., ALMA Observing Tool) for observation details.
For spectral line analysis (e.g., CO or HI), uncomment the spectral cube section and adjust parameters.
The script assumes basic familiarity with CASA tasks and radio astronomy concepts.
For large datasets, consider running CASA on a high-performance computing cluster with mpicasa for parallel processing.

## Contact
For questions, issues, or contributions, contact Mojtaba Raouf at mojtaba.raouf@gmail.com.

Acknowledgments

CASA is developed by the National Radio Astronomy Observatory (NRAO).
This script is inspired by standard radio astronomy data reduction workflows.
