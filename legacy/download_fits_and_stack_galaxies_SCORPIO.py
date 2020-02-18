# -*- coding: utf-8 -*-
# Jose Benavides (jose.astroph@gmail.com)
"""
Documentation:
Este programa recibe información de RA y Dec de un par de galaxias
interactuantes (o cercanas) u otros objetos individual y descarga los datos
.fits correspondientes de datos del Sloan Digital sky Survey (SDSS) de una
lista generada a partir de los filtros u, g, r, i, z. calcula la distancia
en Mpc y la separación (solamente para el caso de pares).
"""
# import pandas as pd
import numpy as np
from astropy import units as apu
from astropy.coordinates import SkyCoord
from astroquery.skyview import SkyView
from retrying import retry
# import matplotlib.image as img
import logging
import os

# logger = logging.getLogger(__scorpio__)
# logger = logging.getLogger(__name__)

"""
1237662225682006144 196.62870697000000 39.844405905521903 0.10972299000000001
1237662225682006156 196.63408457000000 39.849490595521900 0.10918017000000001
AA = [196.63408457000000, 39.849490595521900, 0.10918017000000001]
BB = [196.62870697000000, 39.844405905521903, 0.10972299000000001]

1237651250409767026 126.38991429000001 47.305200665521902 0.12554201000000001
1237651250409767016 126.39162693999999 47.296980665521900 0.12573827000000001
AA = [126.39162693999999, 47.296980665521900, 0.12573827000000001]
BB = [126.38991429000001, 47.305200665521902, 0.12554201000000001]

import scorpio_indv_img as scorpio
scorpio.indv_pair(AA,BB, plx=500)
"""


def stack_pair(glx1, glx2, plx=1000, SURVEY=None, filters=None):
    """
    generate a individual image from list (RA, DEC, z) data of galaxy pair
    for defaul pixles: plx = 1000 and g, i filters
    """
    glx_array = np.array([glx1, glx2])
    df = glx_array
    plx = plx
    # SURVEY = 'SDSS'
    SURVEY = 'SDSS' if SURVEY is None else SURVEY
    filters = ["g", "i"] if filters is None else filters
    optional_filters = ["u", "g", "r", "i", "z"]

    for ff in filters:
        if ff not in optional_filters:
            raise ValueError

    if os.path.exists("missing_filters.dat"):
        os.remove('missing_filters.dat')

    # logger = logging.getLogger('registers_log')
    logger = logging.getLogger('scorpio')
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler('missing_filters.dat')
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)
    fmtLog = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    formatter = logging.Formatter(fmtLog)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # ### function to extract coordinates in Aladin format
    def aladin_coords(pos):
        """function to extract coordinates in Aladin format"""
        print("Coordinates: %2d:%2d:%3.1f %2d:%2d:%3.1f" % (
            int(pos.ra.hms[0]), int(pos.ra.hms[1]), pos.ra.hms[2],
            int(pos.dec.dms[0]), pos.dec.dms[1], pos.dec.dms[2]))
        print(pos.ra.to_string() + " " + pos.dec.to_string())

    # ### Function for downloading images:
    @retry(stop_max_attempt_number=4)
    def download_data(pos):
        path = SkyView.get_images(
            position=pos, survey=SURVEY+str(filters[ff]),
            radius=2*apu.arcmin, pixels=(plx, plx),
            coordinates='J2000', show_progress=True
            )
        return path

    # ############ loop for download images data fits: ################

    # N=len(df)
    missing = []
    stack_glx1 = np.zeros(shape=(plx, plx))
    stack_glx2 = np.zeros(shape=(plx, plx))
    g1g2 = [stack_glx1, stack_glx2]

    for ff in range(len(filters)):
        for ii in range(len(g1g2)):
            pos = SkyCoord(ra=df[ii, 0]*apu.degree, dec=df[ii, 1]*apu.degree)
            aladin_coords(pos)
            try:
                stamp = download_data(pos)
                # image_data = stamp[0][0].data
                # print(stamp[0][0].data)
                g1g2[ii] += stamp[0][0].data

            except Exception:
                missing.append((ii, filters[ff], df[ii, 1], df[ii, 2]))
                # print("Data exist in the directory or",
                # "No data ii:"+str(ii)+" of filter: "+str(filters[ff]))
                # logging.warning("Data exist in the directory or",
                # "No data ii:"+str(ii)+"of filter:"+str(filters[ff]))
                msg1 = "No data filter: "
                logger.info(msg1+str(filters[ff])+" in gal:"+str(ii))
                pass

            print("complete filter:", filters[ff], "of image N:", ii)

    missing = np.array(missing)
    HM = "N, filter,    RA,      DEC"
    np.savetxt("missing_fits.txt", missing, header=HM, delimiter=" ", fmt="%s")

    return g1g2

# END
