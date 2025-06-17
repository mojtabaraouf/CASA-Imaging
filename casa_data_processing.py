
# Title: Radio Astronomy Data Processing Script
# Author: Mojtaba Raouf
# Contact: mojtaba.raouf@gmail.com
# Description: This Python script processes raw radio astronomy data (assumed to be a Measurement Set)
# using CASA (Common Astronomy Software Applications). It performs key signal processing tasks:
# importing data, inspecting for quality, flagging bad data, calibrating instrumental effects,
# imaging the sky, and analyzing results. The script is designed for interferometric data
# (e.g., from ALMA or VLA) and uses a dummy dataset ('dummy_data.ms') for demonstration.
# For real data, replace the dummy MS with an actual Measurement Set and adjust parameters
# (e.g., field names, imaging parameters) based on the observation details.
# Requirements: CASA installed (available from NRAO), Python environment with CASA modules.
# Outputs: Observation summary, diagnostic plots, calibration tables, cleaned image (FITS),
# and image statistics.
# Usage: Run in a CASA environment (e.g., casapy). Update paths and parameters as needed.
# For any details or issues, contact Mojtaba Raouf at mojtaba.raouf@gmail.com.
# Import required CASA modules and standard Python libraries
import os
import glob
from casatasks import listobs, flagdata, bandpass, gaincal, applycal, tclean, imstat, exportfits

# Define paths for dummy data and outputs (replace with actual paths for real data)
# For this example, we assume a dummy Measurement Set file 'dummy_data.ms'
data_dir = './data/'
ms_file = data_dir + 'dummy_data.ms'
output_dir = './outputs/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# ---------------- Section 1: Data Import and Inspection ----------------
# Purpose: Import the raw data (Measurement Set) and inspect its contents to understand
# the observation setup, such as antennas, spectral windows, and sources.
# For dummy data, we assume the MS file exists; in practice, you would use importasdm
# to convert raw telescope data (e.g., ALMA ASDM format) to a Measurement Set.

# Check if the dummy MS file exists (replace with actual import for real data)
if not os.path.exists(ms_file):
    print(f"Warning: {ms_file} not found. Assuming dummy data is available.")
    # For real data, uncomment the following to import ASDM:
    # from casatasks import importasdm
    # importasdm(asdm='raw_data.asdm', vis=ms_file, verbose=True)

# List observation details (antennas, sources, spectral windows, etc.)
listobs(vis=ms_file, listfile=output_dir + 'obs_summary.txt')
print("Observation summary saved to obs_summary.txt")

# Visualize amplitude vs. time to inspect data quality
# PlotMS is typically used interactively, but here we save a plot for reference
from casatools import plotms
plotms(vis=ms_file, xaxis='time', yaxis='amp', plotfile=output_dir + 'amp_vs_time.png', overwrite=True)
print("Amplitude vs. time plot saved to amp_vs_time.png")

# ---------------- Section 2: Data Flagging ----------------
# Purpose: Identify and remove bad data (e.g., RFI, antenna issues) to improve data quality.
# For dummy data, we apply basic flagging rules (e.g., shadow, high amplitudes).
# In practice, inspect plotms output to customize flagging.

# Flag data affected by antenna shadowing (common in interferometry)
flagdata(vis=ms_file, mode='shadow', tolerance=0.0, flagbackup=True)

# Flag high-amplitude outliers (e.g., RFI spikes)
# Adjust threshold based on plotms inspection (here, dummy threshold of 10 Jy)
flagdata(vis=ms_file, mode='clip', clipminmax=[0, 10.0], flagbackup=True)

# Flag specific bad antennas or time ranges (example for dummy data)
# Replace with actual antenna names or times after inspection
flagdata(vis=ms_file, mode='manual', antenna='ANT1', flagbackup=True)

print("Data flagging completed. Backup flags saved.")

# ---------------- Section 3: Calibration ----------------
# Purpose: Correct instrumental and environmental effects (e.g., antenna gains, bandpass).
# For dummy data, we assume calibrator sources are known (e.g., a bandpass calibrator and phase calibrator).
# Replace source names with actual ones from listobs output.

# Bandpass calibration: Correct frequency-dependent effects using a bright calibrator
# Assume 'CAL_BP' is the bandpass calibrator source
bandpass(vis=ms_file, caltable=output_dir + 'bandpass.cal', field='CAL_BP', refant='ANT0',
         solint='inf', combine='scan', minsnr=3.0)
print("Bandpass calibration table generated: bandpass.cal")

# Gain calibration: Correct time-dependent amplitude and phase variations
# Assume 'CAL_PHASE' is the phase calibrator
gaincal(vis=ms_file, caltable=output_dir + 'gain.cal', field='CAL_PHASE', refant='ANT0',
        solint='60s', minsnr=3.0, calmode='ap')
print("Gain calibration table generated: gain.cal")

# Apply calibrations to the science target (assume 'TARGET' is the science source)
applycal(vis=ms_file, field='TARGET', gaintable=[output_dir + 'bandpass.cal', output_dir + 'gain.cal'],
         interp='linear', calwt=False)
print("Calibrations applied to target data.")

# ---------------- Section 4: Imaging ----------------
# Purpose: Create a sky image from calibrated data using the CLEAN algorithm.
# For dummy data, we assume a single field and spectral window.
# Adjust parameters (e.g., cell size, image size) based on actual data properties.

# Create a cleaned image using tclean
tclean(vis=ms_file, imagename=output_dir + 'target_image', field='TARGET', spw='0',
       imsize=[256, 256], cell='1arcsec', niter=1000, threshold='0.1mJy',
       deconvolver='hogbom', weighting='briggs', robust=0.5)
print("Cleaned image generated: target_image")

# Export the image to a FITS file for further analysis
exportfits(imagename=output_dir + 'target_image.image', fitsimage=output_dir + 'target_image.fits')
print("Image exported to target_image.fits")

# ---------------- Section 5: Analysis ----------------
# Purpose: Extract basic statistics from the image (e.g., peak flux, RMS noise).
# For dummy data, we compute image statistics and save them.

# Calculate image statistics (e.g., peak, RMS)
stats = imstat(imagename=output_dir + 'target_image.image')
with open(output_dir + 'image_stats.txt', 'w') as f:
    f.write(f"Peak Flux: {stats['max'][0]:.3f} Jy/beam\n")
    f.write(f"RMS Noise: {stats['rms'][0]:.3f} Jy/beam\n")
print("Image statistics saved to image_stats.txt")

# Optional: Spectral line analysis (if data includes spectral cube)
# For dummy data, assume a spectral cube was created (uncomment for real data)
# tclean(vis=ms_file, imagename=output_dir + 'target_cube', field='TARGET', spw='0',
#        imsize=[256, 256], cell='1arcsec', niter=1000, threshold='0.1mJy',
#        specmode='cube', outframe='LSRK')
# print("Spectral cube generated: target_cube")

# Example: Extract a spectrum at a specific position (replace coordinates with actual ones)
# from casatasks import imval
# spectrum = imval(imagename=output_dir + 'target_cube.image', region='circle [[128pix, 128pix], 5pix]')
# print("Spectrum extracted at central position.")

# ---------------- Section 6: Cleanup ----------------
# Purpose: Organize outputs and clean up temporary files.
# For dummy data, we keep all outputs; in practice, remove unnecessary files.

# List all generated files
print("Generated files:")
for f in glob.glob(output_dir + '*'):
    print(f" - {f}")

# Optional: Remove temporary calibration tables if no longer needed
# os.remove(output_dir + 'bandpass.cal')
# os.remove(output_dir + 'gain.cal')

print("Data processing completed. Outputs saved in", output_dir)