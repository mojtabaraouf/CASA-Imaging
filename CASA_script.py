import numpy as np
import os
import astropy.units as u

# Configuration
pixel = 1024
snaps = ['030', '040', '050', '060', '070', '080']
incl = 41

# Loop through each snapshot to import FITS files
for snap in snaps:
    # Define the FITS files to import
    FIT_IM = [
        f'image_s{snap}_line1_CI_i{incl}_pix{pixel}_Jypixel.fits',
        f'image_s{snap}_line1_CO_i{incl}_pix{pixel}_Jypixel.fits',
    ]

    # Import each FITS file into CASA image format
    for fitsim in FIT_IM:
        importfits(fitsimage=fitsim, imagename=fitsim.split('.')[0] + '.image')

    # Get beam size from the header of the CI line image
    ia.open(f'image_s{snap}_line1_CI_i{incl}_pix{pixel}_Jypixel.image')
    ia.close()

    # Open CO line image to get beam size
    ia.open(f'image_s{snap}_line1_CO_i{incl}_pix{pixel}_Jypixel.image')
    ia.close()

    # List of images for moment generation
    IMG_IM = [
        f'image_s{snap}_line1_CI_i{incl}_pix{pixel}_Jypixel.image',
        f'image_s{snap}_line1_CO_i{incl}_pix{pixel}_Jypixel.image',
    ]

    # Update header information to set units
    for image in IMG_IM:
        imhead(imagename=image, mode='put', hdkey='BUNIT', hdvalue='Jy/pixel')

    # Generate moment maps (0, 1, 2)
    for image in IMG_IM:
        immoments(imagename=image, moments=[0, 1, 2], outfile=image.split('.')[0])

    # Define output extensions and export to FITS format
    ext_casa = ['integrated', 'weighted_coord', 'weighted_dispersion_coord']
    ext_out = ['mom0_Jypixel_kms.fits', 'mom1_kms.fits', 'mom2_kms.fits']

    for ext_c, ext_o in zip(ext_casa, ext_out):
        exportfits(imagename=f'image_s{snap}_line1_CI_i{incl}_pix{pixel}_Jypixel.{ext_c}',
                   fitsimage=f'image_s{snap}_line1_CI_i{incl}_pix{pixel}_{ext_o}')

        exportfits(imagename=f'image_s{snap}_line1_CO_i{incl}_pix{pixel}_Jypixel.{ext_c}',
                   fitsimage=f'image_s{snap}_line1_CO_i{incl}_pix{pixel}_{ext_o}')

# Clean up temporary files
os.system('rm -rf *.last')
os.system('rm -rf *.image *.integrated *.weighted_coord *.weighted_dispersion_coord')