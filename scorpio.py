#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Jose Benavides (jose.astroph@gmail.com)

# This file is part of the
#   SCORPIO Project (https://github.com/josegit88/SCORPIO).
# Copyright (c) 2020, Jose Benavides
# License: MIT
#   Full Text: https://github.com/josegit88/SCORPIO/blob/add-license-2/LICENSE

# ============================================================================
# DOCS
# ============================================================================

"""scorpio."""

# =============================================================================
# IMPORTS
# =============================================================================

import logging
import os
import urllib

import astropy.cosmology as asc
from astropy import units as apu
from astropy import wcs
from astropy.coordinates import SkyCoord

from astroquery.skyview import SkyView

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import numpy as np

from retrying import retry

# ----------------------------------

__version__ = "0.0.1"

logger = logging.getLogger("scorpio")

VALID_FILTERS_SDSS = ["u", "g", "r", "i", "z"]
VALID_FILTERS_2MASS = ["-J", "-H", "-K"]
VALID_FILTERS_WISE = [" 3.4", " 4.6", " 12", " 22"]


# ============================================================================
# CLASSES
# ============================================================================


class NoFilterToStackError(ValueError):
    """Error generated when data of stack galaxies is empty."""

    pass


class NoSurveyInListToStackError(ValueError):
    """Error generated when survey of stack galaxies not in our list."""

    pass


class Image:
    """Main object of the application, receives data from the other functions.

    matrix: applied information from .fits files

    header: header information of some .fits

    dist_physic: estimates the physical distance to the observer in Mpc,
    from a given cosmology

    length_arc: estimates the physical distance between the pair in kpc

    dist_pix: conversion of the physical distance to pixels in the image

    pos1: RA, DEC, Z information of the primary galaxy

    pos2: RA, DEC, Z information of the secondary galaxy

    resolution: integer value of the image resolution in pixels.

    """

    def __init__(
        self,
        matriz=None,
        header=None,
        dist_physic=None,
        dist_pix=None,
        length_arc=None,
        pos1=None,
        pos2=None,
        resolution=None,
    ):
        self.matriz = matriz
        self.header = header
        self.dist_physic = dist_physic
        self.dist_pix = dist_pix
        self.length_arc = length_arc
        self.resolution = resolution
        self.pos1 = pos1
        self.pos2 = pos2

    def plot(
        self,
        ax=None,
        fig=plt.gcf(),
        size=8,
        color_map="inferno",
        dir_images="./individual_images",
        save_Img=False,
        imgName="img1.png",
        **kwargs,
    ):
        r"""Receives data from other functions to generate and export a image.

        ax: Complete information of the properties to generate the image,
            by default it is None.

        dir_images: Destination directory where the images will be exported,
                    by default it is "./individual_images".

        save_Img: Option to save the image, by default it is False.

        imgName: name of the image to be exported, by default it is "img1.png".

        ** kwargs.

        """
        # ------ plot features: ------
        if not os.path.exists(dir_images):
            os.makedirs(dir_images)

        final_imageA = self.matriz
        final_imageA = final_imageA[0]
        plx = self.resolution
        c1 = self.pos1
        c2 = self.pos2
        dis_c1_c2 = self.dist_pix
        s_AB = self.dist_physic

        if ax is None:
            figsize = kwargs.get("figsize", (size, size))
            fig, ax = plt.subplots(figsize=figsize)

        xx = plx / 2.0
        yy = plx / 2.0

        ax.axis([-xx * 0.8, xx * 0.8, -yy * 0.8, yy * 0.8])
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
                if final_imageA[vv, hh] <= 0.0005 * max_value:
                    final_imageA[vv, hh] = 0.0005 * max_value

        max_values_col = []
        min_values_col = []
        for mm in range(len(final_imageA)):
            max_in_column = max(final_imageA[:, mm])
            max_values_col.append(max_in_column)
            min_in_column = min(final_imageA[:, mm])
            min_values_col.append(min_in_column)

        max_value = max(max_values_col)
        min_value = min(min_values_col)

        Norm = mpl.colors.Normalize(vmin=min_value, vmax=max_value)
        Cmap = mpl.cm.ScalarMappable(norm=Norm, cmap=color_map)
        axins1 = inset_axes(
            ax,
            width="5%",
            height="30%",
            loc="lower right",
        )
        cb = fig.colorbar(Cmap, ax=ax, cax=axins1, orientation="vertical")
        cb.set_ticks([])
        ax.imshow(
            final_imageA,
            extent=extent,
            cmap=color_map,
            norm=LogNorm(vmin=min_value, vmax=max_value),
        )
        ax.plot(
            c1[0] - yy,
            -(c1[1] - xx),
            color="cyan",
            marker="o",
            markersize=20,
            mew=2,
            fillstyle="none",
        )
        ax.plot(
            c2[0] - yy,
            -(c2[1] - xx),
            color="cyan",
            marker="o",
            markersize=20,
            mew=2,
            fillstyle="none",
        )

        len_bar = 50.0 * dis_c1_c2 / s_AB
        ax.broken_barh(
            [(-plx / 2.5, len_bar + plx / 7.5)],
            (-plx / 2.5, plx / 7.5),
            facecolors="w",
        )
        ax.hlines(
            y=-plx / 2.72,
            xmin=-plx / 3.0,
            xmax=-plx / 3.0 + len_bar,
            color="k",
            linewidth=3,
        )
        ax.text(
            -plx / 3.0, -plx / 3.0, "50 kpc", fontsize=size * 2.5, color="k"
        )

        if save_Img is True:
            name_image = dir_images + "/" + str(imgName)
            plt.savefig(name_image, bbox_inches="tight", dpi=200)
            plt.close()
        if save_Img is False:
            pass
        # ----------------------------------

        return ax


# -------------

# ============================================================================
# FUNCTIONS
# ============================================================================


@retry(stop_max_attempt_number=4)
def download_data(pos, survey, filters, plx, ff):
    """Download data fits from survey."""
    path = SkyView.get_images(
        position=pos,
        survey=survey + str(filters[ff]),
        radius=2 * apu.arcmin,
        pixels=(plx, plx),
        coordinates="J2000",
        show_progress=True,
    )
    return path


def stack_pair(
    ra1,
    dec1,
    ra2,
    dec2,
    z1=None,
    z2=None,
    resolution=1000,
    survey="SDSS",
    filters=None,
):
    """Generate a individual image from list (RA, DEC, z) data of galaxy pair.

    for defaul pixles: plx = 1000 and g, i filters.

    """
    print(survey)
    plx = resolution

    glx1 = [ra1, dec1, z1]
    glx2 = [ra2, dec2, z2]

    glx_array = np.array([glx1, glx2])

    if survey not in ["SDSS", "2MASS", "WISE"]:
        print("invalid survey")
        raise NoSurveyInListToStackError("Survey not allowed")

    # SDSS:
    if survey == "SDSS":
        filters = ["g", "i"] if filters is None else filters

        for ff in filters:
            if ff not in VALID_FILTERS_SDSS:
                raise NoFilterToStackError(f"{ff} is not a valid filter")

    # 2MASS:
    if survey == "2MASS":
        filters = ["-H", "-K"] if filters is None else filters
        for ff in filters:
            if ff not in VALID_FILTERS_2MASS:
                raise NoFilterToStackError(f"{ff} is not a valid filter")

    # WISE:
    if survey == "WISE":
        filters = [" 4.6", " 12"] if filters is None else filters

        for ff in filters:
            if ff not in VALID_FILTERS_WISE:
                raise NoFilterToStackError(f"{ff} is not a valid filter")

    # ------- download images data fits in diferent filters: --------
    g1g2 = [np.empty(shape=(plx, plx)), np.empty(shape=(plx, plx))]

    for ff in range(len(filters)):
        for ii in range(len(g1g2)):
            pos = SkyCoord(
                ra=glx_array[ii, 0] * apu.degree,
                dec=glx_array[ii, 1] * apu.degree,
            )
            try:
                stamp = download_data(
                    pos=pos, survey=survey, filters=filters, plx=plx, ff=ff
                )
                if ii == 0:
                    stamp_g1 = stamp
            except urllib.error.HTTPError as err:
                if err.code != 404:
                    raise err
                logger.warning(f"No data filter '{filters[ff]}' in gal {ii}")
            else:
                g1g2[ii] += stamp[0][0].data

    if np.all(g1g2[0] == 0):
        raise NoFilterToStackError("Empty array for galaxy1")

    # header = stamp_g1[0][0].header
    header = stamp_g1[0][0].header
    return g1g2, header, plx


def distances(ra1, dec1, ra2, dec2, z1, z2, header, cosmology=asc.Planck15):
    """Receives the RA, DEC, redshift parameters for the two galaxies.

    As well as header information for the primary galaxy and
    cosmology within the options provided by Astropy.

    Calculate the physical and pixel distances of the two galaxies necessary
    for their location in the final image.

    WMAP5        Komatsu et al. 2009            70.2    0.277
    WMAP7        Komatsu et al. 2011            70.4    0.272
    WMAP9        Hinshaw et al. 2013            69.3    0.287
    Planck13     Planck Collab 2013, Paper XVI  67.8    0.307
    Planck15     Planck Collab 2015, Paper XIII 67.7    0.307
    """
    if cosmology not in [
        asc.WMAP5,
        asc.WMAP7,
        asc.WMAP9,
        asc.Planck13,
        asc.Planck15,
    ]:
        raise TypeError("cosmology not allowed")

    glx1 = [ra1, dec1, z1]
    glx2 = [ra2, dec2, z2]

    glx_array = np.array([glx1, glx2])
    z_glx = glx_array[:, 2]
    scale_factor = 1.0 / (1.0 + np.mean(z_glx))

    dist_comv = cosmology.comoving_distance(np.mean(z_glx)).value

    coord_A = SkyCoord(
        ra=glx_array[0, 0] * apu.deg, dec=glx_array[0, 1] * apu.deg
    )
    coord_B = SkyCoord(
        ra=glx_array[1, 0] * apu.deg, dec=glx_array[1, 1] * apu.deg
    )

    theta_rad = coord_A.separation(coord_B).rad
    s_AB = (dist_comv * theta_rad) * 1000.0 * scale_factor

    data_WCS = wcs.WCS(header)

    c1 = data_WCS.wcs_world2pix(glx_array[0, 0], glx_array[0, 1], 0)
    c2 = data_WCS.wcs_world2pix(glx_array[1, 0], glx_array[1, 1], 0)
    dis_c1_c2 = np.sqrt((c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2)

    return s_AB, dis_c1_c2, c1, c2


def gpair(
    ra1,
    dec1,
    ra2,
    dec2,
    z1,
    z2,
    survey="SDSS",
    resolution=1000,
    cosmology=asc.Planck15,
):
    """Receives the RA, DEC, redshift parameters for the two galaxies.

    As well as the resolution in pixels, survey and filters.
    Returns the necessary characteristics to generate the final image.
    """
    img_gp = Image()
    # ---

    g1g2, header, plx = stack_pair(
        ra1,
        dec1,
        ra2,
        dec2,
        z1=z1,
        z2=z2,
        resolution=resolution,
        survey=survey,
    )
    img_gp.matriz = g1g2
    img_gp.header = header
    img_gp.resolution = plx
    # -----

    dist_physic, dist_pix, pos1, pos2 = distances(
        ra1, dec1, ra2, dec2, z1, z2, header, cosmology
    )
    img_gp.dist_physic = dist_physic
    img_gp.dist_pix = dist_pix
    img_gp.pos1 = pos1
    img_gp.pos2 = pos2

    return img_gp


# --------------------

# END
