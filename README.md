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

## Contact

For questions or issues, please contact Mojtaba Raouf(mailto:raouf@strw.leidenuniv.nl).

---

Feel free to modify the README to fit your specific context and needs!
