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

import numpy as np

from astropy import units as apu
from astropy.coordinates import SkyCoord

from astroquery.skyview import SkyView

from retrying import retry


logger = logging.getLogger("scorpio")

VALID_FILTERS = ["u", "g", "r", "i", "z"]


#~ 1237662225682006144 196.62870697000000 39.844405905521903 0.10972299000000001
#~ 1237662225682006156 196.63408457000000 39.849490595521900 0.10918017000000001
#~ AA = [196.63408457000000, 39.849490595521900, 0.10918017000000001]
#~ BB = [196.62870697000000, 39.844405905521903, 0.10972299000000001]

#~ 1237651250409767026 126.38991429000001 47.305200665521902 0.12554201000000001
#~ 1237651250409767016 126.39162693999999 47.296980665521900 0.12573827000000001
#~ AA = [126.39162693999999, 47.296980665521900, 0.12573827000000001]
#~ BB = [126.38991429000001, 47.305200665521902, 0.12554201000000001]

#~ import scorpio_indv_img as scorpio
#~ scorpio.indv_pair(AA,BB, plx=500)


#~ def to_aladin_coords(pos):
    #~ """function to extract coordinates in Aladin format"""
    #~ print("Coordinates: %2d:%2d:%3.1f %2d:%2d:%3.1f" % (
        #~ int(pos.ra.hms[0]), int(pos.ra.hms[1]), pos.ra.hms[2],
        #~ int(pos.dec.dms[0]), pos.dec.dms[1], pos.dec.dms[2]))
    #~ print(pos.ra.to_string() + " " + pos.dec.to_string())

class NoFilterToStackError(RuntimeError):
    """Esfdfsidjsdi cuando

    """
    pass


@retry(stop_max_attempt_number=4)
def download_data(pos, survey, filters, plx):
    """Function for downloading images.

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
    df = glx_array

    filters = ["g", "i"] if filters is None else filters


    for ff in filters:
        if ff not in VALID_FILTERS:
            raise ValueError(f"{ff} is not a valid filter")

    # ############ loop for download images data fits: ################
    g1g2 = [
        np.zeros(shape=(plx, plx)),
        np.zeros(shape=(plx, plx))]

    for ff in range(len(filters)):
        for ii in range(len(g1g2)):
            pos = SkyCoord(ra=df[ii, 0] * apu.degree,
                           dec=df[ii, 1] * apu.degree)
            try:
                stamp = download_data(
                    pos=pos, survey=survey, filters=filters, plx=plx)
            except Exception:
                logger.warning(f"No data filter '{filters[ff]}' in gal {ii}")
            else:
                g1g2[ii] += stamp[0][0].data

    # que pasa si no estakeo un joraca
    if np.all(g1g2 == 0):
        raise NoFilterToStackError("no se encontro ningun filtro para ")
    return g1g2

# END
