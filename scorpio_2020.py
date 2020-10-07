# -*- coding: utf-8 -*-
# Jose Benavides (jose.astroph@gmail.com)

"""Este programa recibe información de RA y Dec de un par de galaxias
interactuantes (o cercanas) u otros objetos individual y descarga los datos
.fits correspondientes de datos del Sloan Digital sky Survey (SDSS) de una
lista generada a partir de los filtros u, g, r, i, z. calcula la distancia
en Mpc y la separación (solamente para el caso de pares).

"""

import logging
import os
import urllib

import numpy as np

from astropy import units as apu
from astropy.coordinates import SkyCoord

from astroquery.skyview import SkyView

from retrying import retry

# add librarys for generated images:
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
plt.close()
# ----------------------------------

__all__ = ["SCORPIO"]
__version__ = "0.0.1"

logger = logging.getLogger("scorpio")

# VALID_FILTERS = ["u", "g", "r", "i", "z"]

VALID_FILTERS_SDSS = ["u", "g", "r", "i", "z"]
VALID_FILTERS_2MASS = ["-J", "-H", "-K"]
VALID_FILTERS_WISE = [" 3.4", " 4.6", " 12", " 22"]


class NoFilterToStackError(RuntimeError):
    """Error generated when data of stack galaxies is empty

    """
    pass


@retry(stop_max_attempt_number=4)
def download_data(pos, survey, filters, plx, ff):
    """Function for download data fits from survey.

    """
    path = SkyView.get_images(
        position=pos, survey=survey + str(filters[ff]),
        radius=2 * apu.arcmin, pixels=(plx, plx),
        coordinates='J2000', show_progress=True)
    return path


def stack_pair(glx1, glx2, plx=1000, survey='SDSS', filters=None):
    """Generate a individual image from list (RA, DEC, z) data of galaxy pair
    for defaul pixles: plx = 1000 and g, i filters.

    """
    glx_array = np.array([glx1, glx2])

    # SDSS:
    if survey == 'SDSS':
        filters = ["g", "i"] if filters is None else filters

        for ff in filters:
            if ff not in VALID_FILTERS_SDSS:
                raise ValueError(f"{ff} is not a valid filter")

    # 2MASS:
    if survey == '2MASS':
        filters = ["-H", "-K"] if filters is None else filters
        for ff in filters:
            if ff not in VALID_FILTERS_2MASS:
                raise ValueError(f"{ff} is not a valid filter")

    # WISE:
    if survey == 'WISE':
        filters = [" 4.6", " 12"] if filters is None else filters

        for ff in filters:
            if ff not in VALID_FILTERS_WISE:
                raise ValueError(f"{ff} is not a valid filter")

    # ------- download images data fits in diferent filters: --------
    g1g2 = [
        np.empty(shape=(plx, plx)),
        np.empty(shape=(plx, plx))]

    for ff in range(len(filters)):
        for ii in range(len(g1g2)):
            pos = SkyCoord(ra=glx_array[ii, 0] * apu.degree,
                           dec=glx_array[ii, 1] * apu.degree)
            try:
                stamp = download_data(
                    pos=pos, survey=survey,
                    filters=filters, plx=plx, ff=ff)
            except urllib.error.HTTPError as err:
                if err.code != 404:
                    raise err
                logger.warning(f"No data filter '{filters[ff]}' in gal {ii}")
            else:
                g1g2[ii] += stamp[0][0].data

    # que pasa si no estakeo un joraca
    if np.all(g1g2[0] == 0):
        raise NoFilterToStackError("Empty array for galaxy1")
    if np.all(g1g2[1] == 0):
        raise NoFilterToStackError("Empty array for galaxy2")
    return g1g2

# generador de imagen:
df = np.array([gal1, gal2]) # array con la info de las dos galaxias

# ====== 1) estimamos la distancia media del observador al par: ========
H0 = 73.52  # constante de Hubble en km/s / Mpc
c_luz = 3.e5  # velocidad de la luz en km/s
D = (c_luz/H0)*np.mean(df[:, 2])  # distance estimated between galaxies [Mpc]
planck = asc.Planck15
dist_comv = planck.comoving_distance(np.mean(df[:, 2])).value
# ======================================================================

# =========== 2) Estima de distancia entre el par: =====================
coord_A = SkyCoord(ra=df[0, 0]*apu.deg, dec=df[0, 1]*apu.deg)
coord_B = SkyCoord(ra=df[1, 0]*apu.deg, dec=df[1, 1]*apu.deg)

theta_rad = coord_A.separation(coord_B).rad
S_AB = (dist_comv*theta_rad)*1000.
# ======================================================================

# ========= 3) Transform from sizes (coords) at pixels: ================
for ff in filters:
    img_test = base_fits+SURVEY+"_image_0_filter_"+str(ff)+".fits" # !!!
    if not os.path.isfile(img_test):  # if not exist the file
        continue
    elif os.path.isfile(img_test):  # if exist the file
        data_imagen = fits.open(img_test)
        break

data_WCS = wcs.WCS(data_imagen[0].header)
# data_imagen[0].header se podria cambiar por stamp[0][0].header

c1 = data_WCS.wcs_world2pix(df[0, 0], df[0, 1], 0)
c2 = data_WCS.wcs_world2pix(df[1, 0], df[1, 1], 0)
dis_c1_c2 = np.sqrt((c1[0]-c2[0])**2+(c1[1]-c2[1])**2)  # pixel-pitagorazed
# =======================================================================

# ============= 4) construccion del plot : ================
dir_images = './individual_images'
erase_options = ["y", "n"]

data_stack = stack_pair(gal1,gal2, plx=100) # esto se debe modificar!!!
final_imageA = data_stack[0]

f, ax = plt.subplots(figsize=(8, 8))

xx = plx/2.
yy = plx/2.

ax.axis([-xx*0.8, xx*0.8, -yy*0.8, yy*0.8])
ax.xaxis.set_major_locator(ticker.NullLocator())
ax.yaxis.set_major_locator(ticker.NullLocator())

extent = [-xx, xx, -yy, yy]

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

save_Img = input("\nYou wish save this image? [y/n]: ")
while save_Img not in erase_options:
    save_Img = input("\nYou wish save this image? [y/n]: ")

if save_Img == "y":
    imgName = input("please input Name for the image: ")
    name_image = dir_ImStack+"/"+str(imgName)
    plt.savefig(name_image, bbox_inches='tight', dpi=200)
    plt.close()
elif save_Img == "n":
    pass

# ======================================================================

# --------------------

# gal1 = [126.39162693999999, 47.296980665521900, 0.12573827000000001]
# gal2 = [126.38991429000001, 47.305200665521902, 0.12554201000000001]

# data_stack2 = stack_pair(gal1,gal2, plx=500)

# END
