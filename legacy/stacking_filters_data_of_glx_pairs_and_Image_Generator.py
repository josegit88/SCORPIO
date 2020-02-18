# -*- coding: utf-8 -*-
# Jose Benavides (jose.astroph@gmail.com)
"""
Documentation:
Este programa recibe los datos .fits de un conjunto de pares de galaxias en
los filtros u, g, r, i, apila los datos y genera una sola imagen incluyendo
una longitud de escala.
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as apu
from astropy.io import fits
from astropy import wcs  # new import
from astropy.coordinates import SkyCoord
import astropy.cosmology as asc
from matplotlib.colors import LogNorm
import matplotlib as mpl
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
plt.close()


# base = "images/images_fits/"
base = "/media/jbenavides/HDD4T/data_images_PISCIS/images/images_fits/"
dir_images_stacked = './images_pairs_stacked'
if not os.path.exists(dir_images_stacked):
    os.popen('mkdir -p ' + dir_images_stacked)

#data_pares = np.genfromtxt("small_sample_data_galaxy_pairs.dat")
data_pares = np.genfromtxt("data_pares_galaxias.dat")

N_pair_generation = 2
plx = 1500

list_galA = range(1, (2*N_pair_generation) + 1, 2)
list_galB = range(0, (2*N_pair_generation), 2)

H0 = 73.52  # constante de Hubble
c_luz = 3.e5  # velocidad de la luz en km/s
D = (c_luz/H0)*data_pares[:, 3]  # distance estimated between galaxies

planck = asc.Planck15
dist_comv = planck.comoving_distance(data_pares[:, 3]).value

imagenes_incompletas = []

for pair in range(N_pair_generation):
    print("pair", pair+1)
    # =============== Estima de distancia entre el par: ===================
    galA = list_galA[pair]
    galB = list_galB[pair]

    # coordenadas de RA y DEC de cada par:
    data_pares[galA, 1], data_pares[galA, 2]
    data_pares[galB, 1], data_pares[galB, 2]

    coord_A = SkyCoord(ra=data_pares[galA, 1]*apu.deg,
                       dec=data_pares[galA, 2]*apu.deg)
    coord_B = SkyCoord(ra=data_pares[galB, 1]*apu.deg,
                       dec=data_pares[galB, 2]*apu.deg)

    theta_rad = coord_A.separation(coord_B).rad
    S_AB = (np.mean((dist_comv[galA], dist_comv[galB]))*theta_rad)*1000.

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

    filtros = ["u", "g", "r", "i"]
    for filtro in filtros:
        # ============= reading and files and stacking: =================
        try:
            ppA = list_galA[pair]
            img_A = base+"SDSS_image_"+str(ppA)+"_filter_"+str(filtro)+".fits"
            imageA_concat.append(fits.getdata(img_A))
        except Exception:
            print("falta el filtro:", filtro, " en la imagen Numero:",
                  list_galA[pair], " del par:", pair+1)
            imagenes_incompletas.append((list_galA[pair], pair+1, filtro))

        try:
            ppB = list_galB[pair]
            img_B = base+"SDSS_image_"+str(ppB)+"_filter_"+str(filtro)+".fits"
            imageB_concat.append(fits.getdata(img_B))
        except Exception:
            print("falta el filtro:", filtro, " en la imagen Numero:",
                  list_galB[pair], " del par:", pair+1)
            imagenes_incompletas.append((list_galB[pair], pair+1, filtro))
        # ===============================================================

    # ===================== sizes at pixels: ========================
    for filtro in filtros:
        ppA = list_galA[pair]
        img_test = base+"SDSS_image_"+str(ppA)+"_filter_"+str(filtro)+".fits"
        # if os.path.isfile(img_test) == False:  # if not exist the file
        if not os.path.isfile(img_test):  # if not exist the file
            continue
        elif os.path.isfile(img_test):  # if exist the file
            data_imagen = fits.open(img_test)
            break

    primer_trozo_imgen = data_imagen[0]
    objeto_WCS = wcs.WCS(primer_trozo_imgen.header)

    c1 = objeto_WCS.wcs_world2pix(data_pares[galA, 1], data_pares[galA, 2], 0)
    c2 = objeto_WCS.wcs_world2pix(data_pares[galB, 1], data_pares[galB, 2], 0)
    dis_c1_c2 = np.sqrt((c1[0]-c2[0])**2+(c1[1]-c2[1])**2)  # pixel-pitagorazed
    # =======================================================================

    # ============== generating the final combined image: =================
    if imageA_concat == []:
        continue
    final_imageA = np.zeros(shape=imageA_concat[0].shape)
    for image in imageA_concat:
        final_imageA += image
    final_imageB = np.zeros(shape=imageB_concat[0].shape)
    for image in imageB_concat:
        final_imageB += image
    # ***********************************************************************

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
    ax.broken_barh([(-plx/2.5, len_bar+plx/7.5)], (-plx/2.5, plx/7.5),
                   facecolors='w')
    ax.hlines(y=-plx/2.72, xmin=-plx/3.0, xmax=-plx/3.0+len_bar, color="k",
              linewidth=3)
    ax.text(-plx/3.0, -plx/3.0, "50 kpc", fontsize=20, color="k")

    name_image = "images_pairs_stacked/galaxy_pair_Number_"+str(pair+1)+".png"
    plt.savefig(name_image, bbox_inches='tight', dpi=200)
    # plt.show()
    plt.close()

imagenes_incompletas = np.array(imagenes_incompletas)
np.savetxt("imagenes_incompletas_con_header.txt", imagenes_incompletas,
           delimiter=" ", header="NrGal NrPair Filter", comments='', fmt="%s")
# np.save("imagenes_incompletas.npy",imagenes_incompletas)
# ---------------------------------
print("-------------------\nThe program has ended successfully")
