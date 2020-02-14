# -*- coding: utf-8 -*-
# Jose Benavides (jose.astroph@gmail.com)
"""
Documentation:
Este programa recibe información de RA y Dec de un conjunto de galaxias
u otros objetos y descarga los datos .fits correspondientes de datos del Sloan
Digital sky Survey (SDSS) e imágenes en formato png en los filtros u, g, r, i.
"""
# import pandas as pd
import numpy as np
from astropy import units as apu
from astropy.coordinates import SkyCoord
from astroquery.skyview import SkyView
from retrying import retry
import matplotlib.image as img
import os


# ### function to extract coordinates in Aladin format
def aladin_coords(pos):
    """function to extract coordinates in Aladin format"""
    print("Coordinates: %2d:%2d:%3.1f %2d:%2d:%3.1f" % (
        int(pos.ra.hms[0]), int(pos.ra.hms[1]), pos.ra.hms[2],
        int(pos.dec.dms[0]), pos.dec.dms[1], pos.dec.dms[2]))
    print(pos.ra.to_string() + " " + pos.dec.to_string())


# ### Load and filter galaxy catalogue:
data_pares = np.genfromtxt('small_sample_data_galaxy_pairs.dat')
df = data_pares  # with some mask

# ### Function for downloading images:
@retry(stop_max_attempt_number=4)
def download_dss(pos):
    path = SkyView.get_images(
        position=pos, survey='SDSS'+str(filters[ff]),
        radius=2*apu.arcmin, pixels=(1500, 1500),
        coordinates='J2000', show_progress=True
        )
    return path


# ### Compute positions, query files and download images (FITS)
base_dir = './'

dir_images = './images'

if not os.path.exists(dir_images):
    os.popen('mkdir -p ' + dir_images)

dir_images_fits = './images/images_fits'
if not os.path.exists(dir_images_fits):
    os.popen('mkdir -p ' + dir_images_fits)

dir_images_png = './images/images_png'
if not os.path.exists(dir_images_png):
    os.popen('mkdir -p ' + dir_images_png)

targets_dir_fits = 'images/images_fits'

targets_dir_fits = 'images/images_fits'
targets_dir_fits = os.path.join(base_dir, targets_dir_fits)

targets_dir_png = 'images/images_png'
targets_dir_png = os.path.join(base_dir, targets_dir_png)


# ############ loop for download images data fits: ################

# N=len(df)

N = 20
filters = ["u", "g", "r", "i"]
missing = []

for ff in range(len(filters)):
    for ii in range(N):
        pos = SkyCoord(ra=df[ii, 1]*apu.degree, dec=df[ii, 2]*apu.degree)
        aladin_coords(pos)
        try:
            stamp = download_dss(pos)
            fits_file = os.path.join(
                targets_dir_fits,
                'SDSS_image_'+str(ii)+'_filter_'+str(filters[ff])+'.fits')
            stamp[0].writeto(fits_file)
            image_data = stamp[0][0].data
            m = image_data.copy()
            m[m < 1.e-5] = 1.e-5
            m = np.log(m)
            png_file = os.path.join(
                targets_dir_png,
                'SDSS_image_'+str(ii)+'_filter_'+str(filters[ff])+'.png')
            img.imsave(png_file, m)
        except Exception:
            missing.append((ii, filters[ff], df[ii, 1], df[ii, 2]))
            print("No image data ii:"+str(ii)+" of filter: "+str(ff))
            pass

        print("complete the filter:", filters[ff], "of image N:", ii)

missing = np.array(missing)
HM = "N, filter,    RA,      DEC"
np.savetxt("missing_fits.txt", missing, header=HM, delimiter=" ", fmt="%s")

# ---------------------------------
print("-------------------\nThe program has ended successfully")
