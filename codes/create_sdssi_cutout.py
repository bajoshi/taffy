from astropy.io import fits
import os

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_dir = home + "/Desktop/ipac/taffy/"

# read in original image and its header
sdss_i = fits.open(taffy_dir + 'SDSS/taffyi_sdss.fits')
hdr = sdss_i[0].header

# define cutout regions
# eyeballed from ds9
"""
 when giving these numbers make sure that the 
 aspect ratio for the original image is preserved.
 i.e. make sure that delta_ra / delta_dec (given in pixels 
 here) is the same for both the new image and the original
"""
dec_pix_low = 280
dec_pix_high = 660
ra_pix_low = 323
ra_pix_high = 703

cutout = sdss_i[0].data[dec_pix_low:dec_pix_high, ra_pix_low:ra_pix_high]
# first dimension is dec i.e. rows in python
# second dimension is ra i.e. cols in python

# replace the original crpix values with the new values
# the rest of the WCS does not change
crpix1 = hdr['CRPIX1'] - ra_pix_low
crpix2 = hdr['CRPIX2'] - dec_pix_low

hdr['CRPIX1'] = crpix1
hdr['CRPIX2'] = crpix2

# rewrite to new file
hdu = fits.PrimaryHDU(data=cutout, header=hdr)
hdu.writeto(taffy_dir + 'SDSS/sdss_i_cutout.fits', clobber=True)