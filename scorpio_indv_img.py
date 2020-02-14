# -*- coding: utf-8 -*-
# Jose Benavides (jose.astroph@gmail.com)
"""
Documentation:
Este programa recibe información de RA y Dec de un par de galaxias
interactuantes (o cercanas) u otros objetos individual y descarga los datos
.fits correspondientes de datos del Sloan Digital sky Survey (SDSS) en los
filtros u, g, r, i. calcula la distancia en Mpc y la separación (solamente
para el caso de pares)
"""
# import pandas as pd
import numpy as np
from astropy import units as apu
from astropy.coordinates import SkyCoord
from astroquery.skyview import SkyView
from retrying import retry
import matplotlib.image as img
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import wcs  # new import
import astropy.cosmology as asc
from matplotlib.colors import LogNorm
import matplotlib as mpl
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from shutil import rmtree
import os
# plt.close()

"""
1237662225682006144 196.62870697000000 39.844405905521903 0.10972299000000001
1237662225682006156 196.63408457000000 39.849490595521900 0.10918017000000001
"""

# ### Load data galaxies:

# glx1 = [196.63408457000000, 39.849490595521900, 0.10918017000000001]
# glx2 = [196.62870697000000, 39.844405905521903, 0.10972299000000001]
# glx_array = np.array([glx1, glx2])
# df = glx_array
# plx = 1500

def indv_pair(glx1,glx2, plx=1500):
    """
    generate a individual image from list (RA, DEC, z) data of galaxy pair
    for defaul pixles: plx = 1500
    """
    glx_array = np.array([glx1, glx2])
    df = glx_array
    # plx = 1500
    plx = plx

    # ### function to extract coordinates in Aladin format
    def aladin_coords(pos):
        """function to extract coordinates in Aladin format"""
        print("Coordinates: %2d:%2d:%3.1f %2d:%2d:%3.1f" % (
            int(pos.ra.hms[0]), int(pos.ra.hms[1]), pos.ra.hms[2],
            int(pos.dec.dms[0]), pos.dec.dms[1], pos.dec.dms[2]))
        print(pos.ra.to_string() + " " + pos.dec.to_string())


    # ### Function for downloading images:
    @retry(stop_max_attempt_number=4)
    def download_dss(pos):
        path = SkyView.get_images(
            position=pos, survey='SDSS'+str(filters[ff]),
            radius=2*apu.arcmin, pixels=(plx, plx),
            coordinates='J2000', show_progress=True
            )
        return path


    # ### Compute positions, query files and download images (FITS)
    base_dir = './'

    dir_images = './individual_images'

    erase_options = ["y", "n"]
    if os.path.exists(dir_images):
        print()
        print("dir /individual_images do exist,"
              "do you wish erase the dir?[y/n]")
        select_dir = input()
        while select_dir not in erase_options:
            print("\ndir /individual_images do exist,"
                  "do you wish erase the dir?[y/n]")
            select_dir = input()
        if select_dir == "y":
            rmtree(dir_images)
        elif select_dir == "n":
            pass

    if not os.path.exists(dir_images):
        os.popen('mkdir -p ' + dir_images)

    dir_images_fits = './individual_images/images_fits'
    if not os.path.exists(dir_images_fits):
        os.popen('mkdir -p ' + dir_images_fits)

    dir_images_png = './individual_images/images_png'
    if not os.path.exists(dir_images_png):
        os.popen('mkdir -p ' + dir_images_png)

    targets_dir_fits = 'individual_images/images_fits'

    targets_dir_fits = 'individual_images/images_fits'
    targets_dir_fits = os.path.join(base_dir, targets_dir_fits)

    targets_dir_png = 'individual_images/images_png'
    targets_dir_png = os.path.join(base_dir, targets_dir_png)


    # ############ loop for download images data fits: ################

    N = len(df)
    filters = ["u", "g", "r", "i"]
    state_download = []

    for ff in range(len(filters)):
        for ii in range(N):
            pos = SkyCoord(ra=df[ii, 0]*apu.degree, dec=df[ii, 1]*apu.degree)
            aladin_coords(pos)
            RA = float("{0:.4f}".format(df[ii, 0]))
            DEC = float("{0:.4f}".format(df[ii, 1]))
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
                state_download.append((ii, filters[ff], RA, DEC, "ok"))
            except Exception:
                state_download.append((ii, filters[ff], RA, DEC, "--"))
                print("No image data ii:"+str(ii)+" of filter: "+str(ff))
                pass

            print("complete the filter:", filters[ff], "of image N:", ii)

    state_download = np.array(state_download)
    print("\nResume:\n   N filter   RA        DEC   state")
    print(state_download)

    # .................... generation of final image: .....................
    base = "individual_images/images_fits/"
    dir_images_stacked = './individual_pairs_stacked'
    if not os.path.exists(dir_images_stacked):
        os.popen('mkdir -p ' + dir_images_stacked)

    
    H0 = 73.52  # constante de Hubble
    c_luz = 3.e5  # velocidad de la luz en km/s
    D = (c_luz/H0)*np.mean(df[:, 2])  # distance estimated between galaxies

    planck = asc.Planck15
    dist_comv = planck.comoving_distance(np.mean(df[:, 2])).value

    imagenes_incompletas = []

    # for pair in range(N_pair_generation):
    # print("pair", pair+1)
    # =============== Estima de distancia entre el par: ===================

    coord_A = SkyCoord(ra=df[0, 0]*apu.deg, dec=df[0, 1]*apu.deg)
    coord_B = SkyCoord(ra=df[1, 0]*apu.deg, dec=df[1, 1]*apu.deg)

    theta_rad = coord_A.separation(coord_B).rad
    S_AB = (dist_comv*theta_rad)*1000.

    # =====================================================================
    f, ax = plt.subplots(figsize=(8, 8))

    llss = 10
    tMl = 5
    tml = 3

    xx = plx/2.
    yy = plx/2.

    ax.axis([-xx*0.8, xx*0.8, -yy*0.8, yy*0.8])
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())

    extent = [-xx, xx, -yy, yy]

    # ********* reading and stacking filters data in only array: ***********
    imageA_concat = []
    imageB_concat = []

    filters = ["u", "g", "r", "i"]
    for ff in filters:
        # ============= reading and files and stacking: =================
        try:
            img_A = base+"SDSS_image_0_filter_"+str(ff)+".fits"
            imageA_concat.append(fits.getdata(img_A))
        except Exception:
            print("falta el filtro:", ff, " en la imagen Numero 0")
            imagenes_incompletas.append((0, 1, ff))

        try:
            img_B = base+"SDSS_image_1_filter_"+str(ff)+".fits"
            imageB_concat.append(fits.getdata(img_B))
        except Exception:
            print("falta el filtro:", ff, " en la imagen Numero 1")
            imagenes_incompletas.append((1, 1, ff))
        # ===============================================================

    # ===================== sizes at pixels: ========================
    for ff in filters:
        img_test = base+"SDSS_image_0_filter_"+str(ff)+".fits"
        if not os.path.isfile(img_test):  # if not exist the file
            continue
        elif os.path.isfile(img_test):  # if exist the file
            data_imagen = fits.open(img_test)
            break

    primer_trozo_imgen = data_imagen[0]
    objeto_WCS = wcs.WCS(primer_trozo_imgen.header)

    c1 = objeto_WCS.wcs_world2pix(df[0, 0], df[0, 1], 0)
    c2 = objeto_WCS.wcs_world2pix(df[1, 0], df[1, 1], 0)
    dis_c1_c2 = np.sqrt((c1[0]-c2[0])**2+(c1[1]-c2[1])**2)  # pixel-pitagoro
    # =======================================================================

    # ============== generating the final combined image: =================
    final_imageA = np.zeros(shape=imageA_concat[0].shape)
    for image in imageA_concat:
        final_imageA += image
    final_imageB = np.zeros(shape=imageB_concat[0].shape)
    for image in imageB_concat:
        final_imageB += image
    # **********************************************************************

    # ax.imshow(final_imageA, extent=extent, cmap='hsv', norm=LogNorm())
    max_values_col = []
    min_values_col = []
    for mm in range(len(final_imageA)):
        max_in_column = max(final_imageA[:, mm])
        max_values_col.append(max_in_column)
        min_in_column = min(final_imageA[:, mm])
        min_values_col.append(min_in_column)

    max_value = max(max_values_col)
    min_value = min(min_values_col)

    for vv in range(len(final_imageA)):
        for hh in range(len(final_imageA)):
            if final_imageA[vv, hh] <= 0.0005*max_value:
                final_imageA[vv, hh] = 0.0005*max_value

    max_values_col = []
    min_values_col = []
    for mm in range(len(final_imageA)):
        max_in_column = max(final_imageA[:, mm])
        max_values_col.append(max_in_column)
        min_in_column = min(final_imageA[:, mm])
        min_values_col.append(min_in_column)

    max_value = max(max_values_col)
    min_value = min(min_values_col)

    norm = mpl.colors.Normalize(vmin=min_value, vmax=max_value)
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.inferno)
    cmap.set_array([])
    axins1 = inset_axes(ax,
                        width="5%",  # width = 50% of parent_bbox width
                        height="30%",  # height : 5%
                        loc='lower right')
    cb = f.colorbar(cmap, ax=ax, cax=axins1, orientation="vertical")
    cb.set_ticks([])
    ax.imshow(final_imageA, extent=extent, cmap='inferno', norm=LogNorm())
    ax.plot(c1[0]-yy, -(c1[1]-xx), color="cyan", marker="o", markersize=20,
            mew=2, fillstyle="none")
    ax.plot(c2[0]-yy, -(c2[1]-xx), color="cyan", marker="o", markersize=20,
            mew=2, fillstyle="none")

    # lscale bar length:
    len_bar = 50.*dis_c1_c2/S_AB
    ax.broken_barh([(-600, len_bar+200)], (-600, 200), facecolors='w')
    ax.hlines(y=-550, xmin=-500, xmax=-500+len_bar, color="k", linewidth=3)
    ax.text(-500, -500, "50 kpc", fontsize=20, color="k")

    imgNum = input("please input Number for image:")
    dir_base = "individual_pairs_stacked"
    name_image = dir_base+"/galaxy_pair_Number_"+str(imgNum)+".png"
    
    plt.savefig(name_image, bbox_inches='tight', dpi=200)
    plt.close()
    
    # ..........................................................................


    # ---- erase directory options: ------
    print("\nYou wish erase fits and png individual files? [y/n]")
    select = input()

    while select not in erase_options:
        print("\nYou wish erase fits and png individual files? [y/n]")
        select = input()

    if select == "y":
        rmtree(dir_images)
    elif select == "n":
        pass
    # -------------------------------------
