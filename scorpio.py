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

__all__ = ["SCORPIO"]
__version__ = "0.0.1"

logger = logging.getLogger("scorpio")

VALID_FILTERS = ["u", "g", "r", "i", "z"]


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

    filters = ["g", "i"] if filters is None else filters


    for ff in filters:
        if ff not in VALID_FILTERS:
            raise ValueError(f"{ff} is not a valid filter")


    # ------- download images data fits in diferent filters: --------
    g1g2 = [
        np.zeros(shape=(plx, plx)),
        np.zeros(shape=(plx, plx))]

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


#gal1 = [126.39162693999999, 47.296980665521900, 0.12573827000000001]
#gal2 = [126.38991429000001, 47.305200665521902, 0.12554201000000001]

#data_stack2 = stack_pair(gal1,gal2, plx=500)

# END
