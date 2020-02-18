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
data_pares = np.genfromtxt('data_pares_galaxias.dat')
df = data_pares  # with some mask
SURVEY = 'SDSS'

# ### Function for downloading images:
@retry(stop_max_attempt_number=4)
def download_data(pos):
    path = SkyView.get_images(
        position=pos, survey=SURVEY+str(filters[ff]),
        radius=2*apu.arcmin, pixels=(plx, plx),
        coordinates='J2000', show_progress=True
        )
    return path


# ### Compute positions, query files and download images (FITS)
base_dir = './'

dir_images = './images'

if not os.path.exists(dir_images):
    os.makedirs(dir_images)

dir_images_fits = './images/images_fits'
if not os.path.exists(dir_images_fits):
    os.makedirs(dir_images_fits)

dir_images_png = './images/images_png'
if not os.path.exists(dir_images_png):
    os.makedirs(dir_images_png)

targets_dir_fits = 'images/images_fits'

targets_dir_fits = 'images/images_fits'
targets_dir_fits = os.path.join(base_dir, targets_dir_fits)

targets_dir_png = 'images/images_png'
targets_dir_png = os.path.join(base_dir, targets_dir_png)


# ############ loop for download images data fits: ################

# N=len(df)

N = 2
plx = 1000
# filters = ["u", "g", "r", "i"]
filters_options = ["u", "g", "r", "i", "z"]
continue_options = ["y", "n"]
filters = []
add_filter = input("select one filter of the list: u, g, r, i, z:\n")
while add_filter not in filters_options:
    print("\nThe filter is not in the options ",
          "u, g, r, i, z")
    add_filter = input("select one filter of the list: u, g, r, i, z:\n")

filters.append(add_filter)
print(filters)
select_opt = input("Do you want to add another filter?[y/n]:\n")
while select_opt not in continue_options:
    select_opt = input("Do you want to add another filter?[y/n]:\n")

while select_opt == "y":
    add_filter = input("select one filter of the list: u, g, r, i, z:\n")
    while add_filter not in filters_options:
        print("\nThe filter is not in the options",
              "u, g, r, i, z")
    if add_filter in filters:
        print("The filter is in the list")
        pass
    elif add_filter not in filters:
        filters.append(add_filter)
    print(filters)
    select_opt = input("Do you want to add another filter?[y/n]:\n")
    while select_opt not in continue_options:
        select_opt = input("Do you want to add another filter?[y/n]:\n")
        if select_opt == "n":
            break

# filters = []
missing = []

image_array = []

for ff in range(len(filters)):
    for ii in range(N):
        pos = SkyCoord(ra=df[ii, 1]*apu.degree, dec=df[ii, 2]*apu.degree)
        aladin_coords(pos)
        try:
            stamp = download_data(pos)
            base_name = SURVEY+'_image_'+str(ii)
            name_fits = '_filter_'+str(filters[ff])+'.fits'
            fits_file = os.path.join(targets_dir_fits, base_name, name_fits)
            stamp[0].writeto(fits_file)
            image_data = stamp[0][0].data
            m = image_data.copy()
            m[m < 1.e-5] = 1.e-5
            m = np.log(m)
            png_file = os.path.join(
                targets_dir_png,
                SURVEY+'_image_'+str(ii)+'_filter_'+str(filters[ff])+'.png')
            img.imsave(png_file, m)
        except Exception:
            missing.append((ii, filters[ff], df[ii, 1], df[ii, 2]))
            print("Data exist in the directory or",
                  "No image data ii:"+str(ii)+" of filter: "+str(filters[ff]))
            pass

        print("complete the filter:", filters[ff], "of image N:", ii)

missing = np.array(missing)
HM = "N, filter,    RA,      DEC"
np.savetxt("missing_fits.txt", missing, header=HM, delimiter=" ", fmt="%s")

# ---------------------------------
print("-------------------\nThe program has ended successfully")
