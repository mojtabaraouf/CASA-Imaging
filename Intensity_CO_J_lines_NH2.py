# HDGAS II Figure 9
#!/usr/bin/env python
"""
Usage:
    intensity_plot.py

The provided Python code generates a scatter plot of intensity and intensity ratios on a FITS image. Specifically, the FITS image is generated using radiative transfer (RADMC-3D) on the HDGAS model. The code requires Python version 3.8 to run and is intended to be executed on the Harbinger platform using the py3 environment.

This file is part of HDGAS simulation.

Copyright (C) 2023 Mojtaba Raouf (mojtaba.raouf@gmail.com)
All Rights Reserved.
"""

import numpy as np
import h5py
from numpy import mean
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.visualization import astropy_mpl_style
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.colors import LogNorm
import pandas
# import pyfits
import copy
from astropy.table import Table
import glob, os
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
import math as ma
import matplotlib.image as mpimg
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import colors
from astropy.utils.data import get_pkg_data_filename
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve
from astropy.stats import sigma_clipped_stats
# from pafit.fit_kinematic_pa import fit_kinematic_pa
# from plotbin.symmetrize_velfield import symmetrize_velfield
# from plotbin.plot_velfield import plot_velfield
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.ticker import AutoMinorLocator
from scipy import stats
from matplotlib import cm
minorLocator = AutoMinorLocator()
visible_ticks = {"top": True, "right": True}

# Path Data- input output
# -----------------------------------------------------------------------------------------------------------
# output = './plots/Intensity/'
output = './plots/CO-dark/'
AGN_model = 'V5B6'
input_AGN_Neq_CI_CO = '/Users/raouf/Work_space/HDGAS/BH3'+str(AGN_model)+'_CH/CI_CO_Neq/'
input_NoAGN_Neq_CI_CO=  '/Users/raouf/Work_space/HDGAS/BH3V0B0_CH/CI_CO_Neq/'

input_AGN_Neq_CII_CO = '/Users/raouf/Work_space/HDGAS/BH3'+str(AGN_model)+'_CH/CII_CO_Neq/'
input_NoAGN_Neq_CII_CO=  '/Users/raouf/Work_space/HDGAS/BH3V0B0_CH/CII_CO_Neq/'

input_AGN_Neq_CI_CII = '/Users/raouf/Work_space/HDGAS/BH3'+str(AGN_model)+'_CH/CI_CII_Neq/'
input_NoAGN_Neq_CI_CII=  '/Users/raouf/Work_space/HDGAS/BH3V0B0_CH/CI_CII_Neq/'

# CO (2-1)
input_AGN_Neq_CI_CO12 = '/Users/raouf/Work_space/HDGAS/BH3'+str(AGN_model)+'_CH/CI_CO21_Neq/'
input_NoAGN_Neq_CI_CO12=  '/Users/raouf/Work_space/HDGAS/BH3V0B0_CH/CI_CO21_Neq/'

input_AGN_Neq_CII_CO12 = '/Users/raouf/Work_space/HDGAS/BH3'+str(AGN_model)+'_CH/CII_CO21_Neq/'
input_NoAGN_Neq_CII_CO12=  '/Users/raouf/Work_space/HDGAS/BH3V0B0_CH/CII_CO21_Neq/'

input_AGN_Neq_CO23 = '/Users/raouf/Work_space/HDGAS/BH3'+str(AGN_model)+'_CH/CO_J_lines/'
input_NoAGN_Neq_CO23=  '/Users/raouf/Work_space/HDGAS/BH3'+str(AGN_model)+'_CH/CO_J_lines/'

input_AGN_Neq_HCOp_CO = '/Users/raouf/Work_space/HDGAS/BH3'+str(AGN_model)+'_CH/HCOp_CO_Neq/'
input_NoAGN_Neq_HCOp_CO=  '/Users/raouf/Work_space/HDGAS/BH3V0B0_CH/HCOp_CO_Neq/'

input_AGN_Neq_H2O_CO = '/Users/raouf/Work_space/HDGAS/BH3'+str(AGN_model)+'_CH/H2O_HCOp_Neq/'
input_NoAGN_Neq_H2O_CO=  '/Users/raouf/Work_space/HDGAS/BH3V0B0_CH/H2O_HCOp_Neq/'

input_AGN_Neq_H2O_CO = '/Users/raouf/Work_space/HDGAS/BH3'+str(AGN_model)+'_CH/H2O_CO_Neq/'
input_NoAGN_Neq_H2O_CO=  '/Users/raouf/Work_space/HDGAS/BH3V0B0_CH/H2O_CO_Neq/'
input_NH2_AGN = '/Users/raouf/Work_space/HDGAS/Snapshots/BH3'+str(AGN_model)+'_CH/'
input_NH2_NoAGN = '/Users/raouf/Work_space/HDGAS/Snapshots/BH3V0B0_CH/'


fn_NH2_CO = [
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_030.hdf5',
                input_AGN_Neq_CI_CO + "image_s030_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_040.hdf5',
                input_AGN_Neq_CI_CO + "image_s040_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_050.hdf5',
                input_AGN_Neq_CI_CO + "image_s050_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_060.hdf5',
                input_AGN_Neq_CI_CO + "image_s060_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_070.hdf5',
                input_AGN_Neq_CI_CO + "image_s070_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_080.hdf5',
                input_AGN_Neq_CI_CO + "image_s080_line1_CO_i41_pix1024_mom0_K_kms.fits",

                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_030.hdf5',
                input_NoAGN_Neq_CI_CO + "image_s030_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_040.hdf5',
                input_NoAGN_Neq_CI_CO + "image_s040_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_050.hdf5',
                input_NoAGN_Neq_CI_CO + "image_s050_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_060.hdf5',
                input_NoAGN_Neq_CI_CO + "image_s060_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_070.hdf5',
                input_NoAGN_Neq_CI_CO + "image_s070_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_080.hdf5',
                input_NoAGN_Neq_CI_CO + "image_s080_line1_CO_i41_pix1024_mom0_K_kms.fits",
                ]

fn_NH2_CO12 = [
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_030.hdf5',
                input_AGN_Neq_CI_CO12 + "image_s030_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_040.hdf5',
                input_AGN_Neq_CI_CO12 + "image_s040_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_050.hdf5',
                input_AGN_Neq_CI_CO12 + "image_s050_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_060.hdf5',
                input_AGN_Neq_CI_CO12 + "image_s060_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_070.hdf5',
                input_AGN_Neq_CI_CO12 + "image_s070_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_080.hdf5',
                input_AGN_Neq_CI_CO12 + "image_s080_line1_CO_i41_pix1024_mom0_K_kms.fits",

                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_030.hdf5',
                input_NoAGN_Neq_CI_CO12 + "image_s030_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_040.hdf5',
                input_NoAGN_Neq_CI_CO12 + "image_s040_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_050.hdf5',
                input_NoAGN_Neq_CI_CO12 + "image_s050_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_060.hdf5',
                input_NoAGN_Neq_CI_CO12 + "image_s060_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_070.hdf5',
                input_NoAGN_Neq_CI_CO12 + "image_s070_line1_CO_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_080.hdf5',
                input_NoAGN_Neq_CI_CO12 + "image_s080_line1_CO_i41_pix1024_mom0_K_kms.fits",
                ]
fn_NH2_CO23 = [
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_030.hdf5',
                input_AGN_Neq_CO23 + "image_s030_line3_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_040.hdf5',
                input_AGN_Neq_CO23 + "image_s040_line3_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_050.hdf5',
                input_AGN_Neq_CO23 + "image_s050_line3_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_060.hdf5',
                input_AGN_Neq_CO23 + "image_s060_line3_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_070.hdf5',
                input_AGN_Neq_CO23 + "image_s070_line3_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_080.hdf5',
                input_AGN_Neq_CO23 + "image_s080_line3_i41_pix1024_mom0_K_kms.fits",

                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_030.hdf5',
                input_NoAGN_Neq_CO23 + "image_s030_line3_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_040.hdf5',
                input_NoAGN_Neq_CO23 + "image_s040_line3_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_050.hdf5',
                input_NoAGN_Neq_CO23 + "image_s050_line3_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_060.hdf5',
                input_NoAGN_Neq_CO23 + "image_s060_line3_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_070.hdf5',
                input_NoAGN_Neq_CO23 + "image_s070_line3_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_080.hdf5',
                input_NoAGN_Neq_CO23 + "image_s080_line3_i41_pix1024_mom0_K_kms.fits",
                ]

fn_NH2_CO34 = [
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_030.hdf5',
                input_AGN_Neq_CO23 + "image_s030_line4_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_040.hdf5',
                input_AGN_Neq_CO23 + "image_s040_line4_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_050.hdf5',
                input_AGN_Neq_CO23 + "image_s050_line4_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_060.hdf5',
                input_AGN_Neq_CO23 + "image_s060_line4_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_070.hdf5',
                input_AGN_Neq_CO23 + "image_s070_line4_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_080.hdf5',
                input_AGN_Neq_CO23 + "image_s080_line4_i41_pix1024_mom0_K_kms.fits",

                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_030.hdf5',
                input_NoAGN_Neq_CO23 + "image_s030_line4_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_040.hdf5',
                input_NoAGN_Neq_CO23 + "image_s040_line4_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_050.hdf5',
                input_NoAGN_Neq_CO23 + "image_s050_line4_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_060.hdf5',
                input_NoAGN_Neq_CO23 + "image_s060_line4_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_070.hdf5',
                input_NoAGN_Neq_CO23 + "image_s070_line4_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_080.hdf5',
                input_NoAGN_Neq_CO23 + "image_s080_line4_i41_pix1024_mom0_K_kms.fits",
                ]
fn_NH2_CO45 = [
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_030.hdf5',
                input_AGN_Neq_CO23 + "image_s030_line5_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_040.hdf5',
                input_AGN_Neq_CO23 + "image_s040_line5_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_050.hdf5',
                input_AGN_Neq_CO23 + "image_s050_line5_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_060.hdf5',
                input_AGN_Neq_CO23 + "image_s060_line5_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_070.hdf5',
                input_AGN_Neq_CO23 + "image_s070_line5_i41_pix1024_mom0_K_kms.fits",
                input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_080.hdf5',
                input_AGN_Neq_CO23 + "image_s080_line5_i41_pix1024_mom0_K_kms.fits",

                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_030.hdf5',
                input_NoAGN_Neq_CO23 + "image_s030_line5_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_040.hdf5',
                input_NoAGN_Neq_CO23 + "image_s040_line5_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_050.hdf5',
                input_NoAGN_Neq_CO23 + "image_s050_line5_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_060.hdf5',
                input_NoAGN_Neq_CO23 + "image_s060_line5_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_070.hdf5',
                input_NoAGN_Neq_CO23 + "image_s070_line5_i41_pix1024_mom0_K_kms.fits",
                input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_080.hdf5',
                input_NoAGN_Neq_CO23 + "image_s080_line5_i41_pix1024_mom0_K_kms.fits",
                ]

fns_ci = [
    input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_030.hdf5',
    input_AGN_Neq_CI_CII + "image_s030_line1_CI_i41_pix1024_mom0_K_kms.fits",
    input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_040.hdf5',
    input_AGN_Neq_CI_CII + "image_s040_line1_CI_i41_pix1024_mom0_K_kms.fits",
    input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_050.hdf5',
    input_AGN_Neq_CI_CII + "image_s050_line1_CI_i41_pix1024_mom0_K_kms.fits",
    input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_060.hdf5',
    input_AGN_Neq_CI_CII + "image_s060_line1_CI_i41_pix1024_mom0_K_kms.fits",
    input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_070.hdf5',
    input_AGN_Neq_CI_CII + "image_s070_line1_CI_i41_pix1024_mom0_K_kms.fits",
    input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_080.hdf5',
    input_AGN_Neq_CI_CII + "image_s080_line1_CI_i41_pix1024_mom0_K_kms.fits",
    input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_030.hdf5',
    input_NoAGN_Neq_CI_CII + "image_s030_line1_CI_i41_pix1024_mom0_K_kms.fits",
    input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_040.hdf5',
    input_NoAGN_Neq_CI_CII + "image_s040_line1_CI_i41_pix1024_mom0_K_kms.fits",
    input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_050.hdf5',
    input_NoAGN_Neq_CI_CII + "image_s050_line1_CI_i41_pix1024_mom0_K_kms.fits",
    input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_060.hdf5',
    input_NoAGN_Neq_CI_CII + "image_s060_line1_CI_i41_pix1024_mom0_K_kms.fits",
    input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_070.hdf5',
    input_NoAGN_Neq_CI_CII + "image_s070_line1_CI_i41_pix1024_mom0_K_kms.fits",
    input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_080.hdf5',
    input_NoAGN_Neq_CI_CII + "image_s080_line1_CI_i41_pix1024_mom0_K_kms.fits",
      ]
fns_cii = [
    input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_030.hdf5',
    input_AGN_Neq_CI_CII + "image_s030_line1_CII_i41_pix1024_mom0_K_kms.fits",
    input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_040.hdf5',
    input_AGN_Neq_CI_CII + "image_s040_line1_CII_i41_pix1024_mom0_K_kms.fits",
    input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_050.hdf5',
    input_AGN_Neq_CI_CII + "image_s050_line1_CII_i41_pix1024_mom0_K_kms.fits",
    input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_060.hdf5',
    input_AGN_Neq_CI_CII + "image_s060_line1_CII_i41_pix1024_mom0_K_kms.fits",
    input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_070.hdf5',
    input_AGN_Neq_CI_CII + "image_s070_line1_CII_i41_pix1024_mom0_K_kms.fits",
    input_NH2_AGN+ 'data_AMRmaps_H2CO_N1024_080.hdf5',
    input_AGN_Neq_CI_CII + "image_s080_line1_CII_i41_pix1024_mom0_K_kms.fits",
    input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_030.hdf5',
    input_NoAGN_Neq_CI_CII + "image_s030_line1_CII_i41_pix1024_mom0_K_kms.fits",
    input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_040.hdf5',
    input_NoAGN_Neq_CI_CII + "image_s040_line1_CII_i41_pix1024_mom0_K_kms.fits",
    input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_050.hdf5',
    input_NoAGN_Neq_CI_CII + "image_s050_line1_CII_i41_pix1024_mom0_K_kms.fits",
    input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_060.hdf5',
    input_NoAGN_Neq_CI_CII + "image_s060_line1_CII_i41_pix1024_mom0_K_kms.fits",
    input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_070.hdf5',
    input_NoAGN_Neq_CI_CII + "image_s070_line1_CII_i41_pix1024_mom0_K_kms.fits",
    input_NH2_NoAGN+ 'data_AMRmaps_H2CO_N1024_080.hdf5',
    input_NoAGN_Neq_CI_CII + "image_s080_line1_CII_i41_pix1024_mom0_K_kms.fits",
      ]
snaps = [30,30,40,40,50,50,60,60,70,70,80,80,30,30,40,40,50,50,60,60,70,70,80,80]

# Functions
# ---------------------------------------------------------------------------------
# color map
def reverse_colourmap(cmap, name = 'my_cmap_r'):
    reverse = []
    k = []

    for key in cmap._segmentdata:
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:
            data.append((1-t[0],t[2],t[1]))
        reverse.append(sorted(data))

    LinearL = dict(zip(k,reverse))
    my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL)
    return my_cmap_r
# color_map = plt.cm.get_cmap('RdYlBu')
# color_map = plt.cm.get_cmap('jet')
color_map = plt.cm.get_cmap('Blues')
#color_map = plt.cm.get_cmap('hot')

#color_map = plt.cm.get_cmap('rainbow')
reversed_color_map = reverse_colourmap(color_map)
# cmap=reversed_color_map
cmap=color_map
def deg_to_parsec(cdelt2):
    """
    Converts the pixel scale in degrees per pixel to parsecs per pixel.

    Parameters:
    cdelt2 (float): The pixel scale in the y-direction (degrees/pixel).

    Returns:
    float: The pixel scale in parsecs/pixel.
    """
    # Convert degrees to radians
    deg_to_rad = np.pi / 180

    # Convert radians to parsecs
    rad_to_parsec = 1 / 206265

    # Calculate the pixel scale in parsecs/pixel
    parsec_per_pixel = cdelt2 * deg_to_rad * rad_to_parsec

    return parsec_per_pixel

def calculate_radii_in_parsec(x0, y0, map_width_pixels, map_height_pixels, cdelt1, cdelt2, distance_mpc, inclination_deg):
    """
    Calculates the radii in parsecs for each pixel in the map.
    Mojtaba Raouf
    Parameters:
    x0 (int): The x-coordinate of the central pixel.
    y0 (int): The y-coordinate of the central pixel.
    map_width_pixels (int): The width of the map in pixels.
    map_height_pixels (int): The height of the map in pixels.
    cdelt1 (float): The pixel scale in the x-direction (degrees/pixel).
    cdelt2 (float): The pixel scale in the y-direction (degrees/pixel).
    distance_mpc (float): The distance to the object in Megaparsecs.
    inclination_deg (float): The inclination angle in degrees.

    Returns:
    list: A list of radii in parsecs for each pixel in the map.
    """
    # Convert inclination from degrees to radians
    inclination_rad = np.deg2rad(inclination_deg)

    # Initialize list to store radii in parsecs
    radii_in_parsec = []

    # Loop through each pixel in the map
    for i in range(map_width_pixels):
        for j in range(map_height_pixels):
            # Calculate radial distance from the center using pixel scale (cdelt1 and cdelt2)
            r = np.sqrt((i - x0)**2 + (1 / np.cos(inclination_rad) * (j - y0))**2) * cdelt2

            # Convert radius from degrees to radians
            r_rad = np.deg2rad(r)

            # Convert radius from radians to parsecs
            radius_in_parsec = r_rad * distance_mpc * 1000000

            # Append to radii list
            radii_in_parsec.append(radius_in_parsec)

    return radii_in_parsec

#Radial profile
def radial_profile(data, center):
    x, y = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile


# Mask with diffrent radii
def create_circular_mask(h, w, center=None, radius=None):

    if center is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask

def runningmedian(x,y,xlolim=-1.e20,ylolim=-1.e20,bins=10,stat='median'):
        xp = x[(x>xlolim)&(y>ylolim)]
        yp = y[(x>xlolim)&(y>ylolim)]
        if bins < 0:	# bins<0 sets bins such that there are equal numbers per bin
            bin_edges = histedges_equalN(xp,-bins)
            bin_means, bin_edges, binnumber = stats.binned_statistic(xp,yp,bins=bin_edges,statistic=stat)
        else:
            bin_means, bin_edges, binnumber = stats.binned_statistic(xp,yp,bins=bins,statistic=stat)
        bin_cent = 0.5*(bin_edges[1:]+bin_edges[:-1])
        ymed = []
        ymean = []
        ysigma = []
        for i in range(0,len(bin_edges[:-1])):
                xsub = xp[xp>bin_edges[i]]
                ysub = yp[xp>bin_edges[i]]
                ysub = ysub[xsub<bin_edges[i+1]]
                ymed.append(np.median(10**ysub))
                ymean.append(np.mean(10**ysub))
                ysigma.append(np.std(10**ysub))
        if stat=='median': ymean = np.asarray(ymed)
        else: ymean = np.asarray(ymean)
        ysiglo = np.maximum(ymean-ysigma,ymean*0.1)
        ysiglo = np.log10(ymean)-np.log10(ysiglo)
        ysighi = np.log10(ymean+ysigma)-np.log10(ymean)
        ymean = np.log10(ymean)
        #print bin_cent,ymean,ysiglo,ysighi
        #plt.plot(bin_cent,ymed,'ro',ms=12,color='c')
        #plt.plot(bin_cent,ymean,'--',lw=3,color='m')
        #plt.errorbar(bin_cent,ymean,yerr=[ysiglo,ysighi],fmt='ro')
        return bin_cent,ymean,ysiglo,ysighi
#Scatter plot from intensity maps
plt.figure(figsize=(8,6))

def create_scatter_plot(image1, image2, xlabel, ylabel, name, rad50, rad100, xmin1, xmin2, ymin1, ymin2, snap, model, xi=None, xj=None, yi=None, yj=None):
    x1, x2, y1, y2 = xmin1, xmin2, ymin1, ymin2
    res = 512

    # CO & H2 column density
    # -------------------------------------------------------
    # fn_AGN =  '/data2/mojtaba/gizmo-public/output_Disk_0103_AGN/BH3'+str(AGN_model)+'_CH/data_AMRmaps_H2CO_N1024_%03d.hdf5' % snap
    f_agn = h5py.File(image1, 'r')
    df_CO_agn = f_agn['AMRmaps']['NCO'][:]
    df_H2_agn = f_agn['AMRmaps']['NH2'][:]
    UnitLength_in_cm   = 3.085678e21 # / 3e18 for convert unit to pc-2
    X=np.ones(shape=f_agn['AMRmaps']['NCO'].shape,dtype=np.float64)
    X *= UnitLength_in_cm
    df_CO_agn = df_CO_agn * X
    df_H2_agn = df_H2_agn * X
    Map_NH2_50 = df_H2_agn.copy()
    Map_NH2_100 = df_H2_agn.copy()
    Map_NH2_300 = df_H2_agn.copy()
    # print(df_H2_agn.shape())

    hdu_2 = fits.open(image2, memmap=True)

    Map_2 = hdu_2[0].data
    Map_H2_50 = Map_2.copy()
    Map_H2_100 = Map_2.copy()
    Map_H2_300 = Map_2.copy()
    z2, h2, w2 = Map_2.shape
    mask_vel2_50 = create_circular_mask(h2, w2, center=(res, res), radius=rad50)
    mask_vel2_100 = create_circular_mask(h2, w2, center=(res, res), radius=rad100)
    mask_vel2_300 = create_circular_mask(h2, w2, center=(res, res), radius=500)
    Map_H2_50[-1, ~mask_vel2_50] = np.nan
    Map_H2_100[-1, ~mask_vel2_100] = np.nan
    Map_H2_300[-1, ~mask_vel2_300] = np.nan

    Map_NH2_300[~mask_vel2_300] = np.nan

    Map_NH2_100[~mask_vel2_100] = np.nan
    Map_NH2_50[~mask_vel2_50] = np.nan
    time = snap / 10.
    print(time, name)
    # x = np.log10(Map_H1_300.flatten())
    y = np.log10(Map_H2_100.flatten())
    z = np.log10(Map_NH2_100.flatten())

    radii_in_parsec = calculate_radii_in_parsec(x0=512, y0=512, map_width_pixels=1024, map_height_pixels=1024, cdelt1=-1.998317230538E-06, cdelt2=1.998317230538E-06, distance_mpc=14, inclination_deg=41)

    # plt.scatter(z, y-z, c=radii_in_parsec, cmap='jet', marker='.', alpha=0.5, s=2, label='t = %0d Myr' % time)
    # cbar = plt.colorbar(label='r [pc]')

    # bin_cent,ymean,ysiglo,ysighi = runningmedian(z, y-z,xlolim=10,ylolim=-26,bins=20,stat='median')
    # plt.plot(bin_cent,ymean,color = 'red',lw=2,alpha=.99, markersize=1, label='median')
    # Upe = ymean + ysighi
    # Lwo = ymean - ysiglo
    # plt.fill_between(bin_cent, Upe, Lwo, facecolor='grey', alpha=0.4)
    bin_cent,ymean,ysiglo,ysighi = runningmedian(z, y,xlolim=16,ylolim=-2,bins=10,stat='median')

    xx = [0,1,2]
    yy = [0,1,1]
    if model == 'AGN':
        c1 = next(color1)
        plt.plot(bin_cent,ymean,color = c1,lw=3,alpha=.99, markersize=1, label='t = %0d Myr' % time)
        if snap == 80:
            plt.plot(xx,yy,'k-', label = 'AGN')
    elif model == 'NoAGN':
        c2 = next(color2)
        plt.plot(bin_cent,ymean,color = c2,lw=3,ls='--',alpha=.99, markersize=1)
        if snap == 80:
            plt.plot(xx,yy,'k--', label = 'NoAGN')

    if xi is not None:
        obsx = [xi, xj]
        obsy = [yi, yj]
        obseq = [0, 0]
        plt.plot(obsx, obsy, 'k:', lw=2)
        plt.plot(obsx, obseq, 'k:', lw=2)

    # plt.text(0.2, 0.8, model, size=20, transform=plt.gcf().transFigure)
    plt.xlabel(xlabel, size=22)
    plt.ylabel(ylabel, size=22)
    #
    plt.xlim(x1, x2)
    plt.ylim(y1, y2)
    # plt.xscale('log')
    # plt.yscale('log')
    plt.subplots_adjust(wspace=0.1, hspace=0.35, bottom=0.15, right=0.96, left=0.15, top=0.98)

    plt.tick_params(which='both', direction='in', **visible_ticks)
    plt.tick_params(labelsize=20)


# plots
# ------------------------------------------------------------------------------------------------------------
j=0
plt.figure(figsize=(8,6))

color1 = iter(cm.rainbow(np.linspace(0, 1, 6)))
color2 = iter(cm.rainbow(np.linspace(0, 1, 6)))
for i in range(0, 24, 2):
    print(i)
    if i < 12:
        create_scatter_plot(fn_NH2_CO45[i], fn_NH2_CO45[i+1], r'$\log\ (\rm N_{H_2}\ [cm^{-2}])$', r'$\log\ (\rm I_{CO(5-4)} [K\ km\ s^{-1}])$', 'CO12_H2_AGN_all',50,100, 19,24,-2,4,snaps[i],'AGN')
    else:
        create_scatter_plot(fn_NH2_CO45[i], fn_NH2_CO45[i+1], r'$\log\ (\rm N_{H_2}\ [cm^{-2}])$', r'$\log\ (\rm I_{CO(5-4)} [K\ km\ s^{-1}])$', 'CO45_H2_NoAGN_all',50,100, 19,24,-2,4,snaps[i],'NoAGN')


# plt.subplots_adjust(wspace=0.1, hspace=0.35, bottom=0.15, right=0.96, left=0.15, top=0.98)

filename = output + 'Intensity_'+str(AGN_model)+'_CO45_H2.png'
plt.savefig(filename)
plt.close()

color1 = iter(cm.rainbow(np.linspace(0, 1, 6)))
color2 = iter(cm.rainbow(np.linspace(0, 1, 6)))
for i in range(0, 24, 2):
    print(i)
    if i < 12:
        create_scatter_plot(fn_NH2_CO34[i], fn_NH2_CO34[i+1], r'$\log\ (\rm N_{H_2}\ [cm^{-2}])$', r'$\log\ (\rm I_{CO(4-3)} [K\ km\ s^{-1}])$', 'CO12_H2_AGN_all',50,100, 19,24,-2,4,snaps[i],'AGN')
    else:
        create_scatter_plot(fn_NH2_CO34[i], fn_NH2_CO34[i+1], r'$\log\ (\rm N_{H_2}\ [cm^{-2}])$', r'$\log\ (\rm I_{CO(4-3)} [K\ km\ s^{-1}])$', 'CO34_H2_NoAGN_all',50,100, 19,24,-2,4,snaps[i],'NoAGN')
# plt.subplots_adjust(wspace=0.1, hspace=0.35, bottom=0.15, right=0.96, left=0.15, top=0.88)

filename = output + 'Intensity_'+str(AGN_model)+'_CO34_H2.png'
plt.savefig(filename)
plt.close()

plt.figure(figsize=(8,6))
color1 = iter(cm.rainbow(np.linspace(0, 1, 6)))
color2 = iter(cm.rainbow(np.linspace(0, 1, 6)))
for i in range(0, 24, 2):
    print(i)
    if i < 12:
        create_scatter_plot(fn_NH2_CO23[i], fn_NH2_CO23[i+1], r'$\log\ (\rm N_{H_2}\ [cm^{-2}])$', r'$\log\ (\rm I_{CO(3-2)} [K\ km\ s^{-1}])$', 'CO12_H2_AGN_all',50,100, 19,24,-2,4,snaps[i],'AGN')
    else:
        create_scatter_plot(fn_NH2_CO23[i], fn_NH2_CO23[i+1], r'$\log\ (\rm N_{H_2}\ [cm^{-2}])$', r'$\log\ (\rm I_{CO(3-2)} [K\ km\ s^{-1}])$', 'CO23_H2_NoAGN_all',50,100, 19,24,-2,4,snaps[i],'NoAGN')
# plt.subplots_adjust(wspace=0.1, hspace=0.35, bottom=0.15, right=0.96, left=0.15, top=0.88)

filename = output + 'Intensity_'+str(AGN_model)+'_CO23_H2.png'
plt.savefig(filename)
plt.close()

plt.figure(figsize=(8,6))
color1 = iter(cm.rainbow(np.linspace(0, 1, 6)))
color2 = iter(cm.rainbow(np.linspace(0, 1, 6)))
for i in range(0, 24, 2):
    print(i)
    if i < 12:
        create_scatter_plot(fn_NH2_CO12[i], fn_NH2_CO12[i+1], r'$\log\ (\rm N_{H_2}\ [cm^{-2}])$', r'$\log\ (\rm I_{CO(2-1)} [K\ km\ s^{-1}])$', 'CO12_H2_AGN_all',50,100, 19,24,-2,4,snaps[i],'AGN')
    else:
        create_scatter_plot(fn_NH2_CO12[i], fn_NH2_CO12[i+1], r'$\log\ (\rm N_{H_2}\ [cm^{-2}])$', r'$\log\ (\rm I_{CO(2-1)} [K\ km\ s^{-1}])$', 'CO12_H2_NoAGN_all',50,100, 19,24,-2,4,snaps[i],'NoAGN')
# plt.subplots_adjust(wspace=0.1, hspace=0.35, bottom=0.15, right=0.96, left=0.15, top=0.88)

filename = output + 'Intensity_'+str(AGN_model)+'_CO12_H2.png'
plt.savefig(filename)
plt.close()

plt.figure(figsize=(8,6))

color1 = iter(cm.rainbow(np.linspace(0, 1, 6)))
color2 = iter(cm.rainbow(np.linspace(0, 1, 6)))
for i in range(0, 24, 2):
    print(i)
    if i < 12:
        create_scatter_plot(fn_NH2_CO[i], fn_NH2_CO[i+1], r'$\log\ (\rm N_{H_2}\ [cm^{-2}])$', r'$\log\ (\rm I_{CO(1-0)} [K\ km\ s^{-1}])$', 'CO_H2_AGN_all',50,100, 19,24,-2,4,snaps[i],'AGN')
    else:
        create_scatter_plot(fn_NH2_CO[i], fn_NH2_CO[i+1], r'$\log\ (\rm N_{H_2}\ [cm^{-2}])$', r'$\log\ (\rm I_{CO(1-0)} [K\ km\ s^{-1}])$', 'CO_H2_NoAGN_all',50,100, 19,24,-2,4,snaps[i],'NoAGN')
leg = plt.legend(loc='lower left',ncol=2)
for t in leg.get_texts():
    t.set_fontsize(30)
# plt.subplots_adjust(wspace=0.1, hspace=0.35, bottom=0.15, right=0.96, left=0.15, top=0.88)

filename = output + 'Intensity_'+str(AGN_model)+'_CO_H2.png'
plt.savefig(filename)
plt.close()
