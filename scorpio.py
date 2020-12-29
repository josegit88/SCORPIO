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

"""Sky COllector of galaxy Pairs and Image Output (Scorpio).

Is a tool to quick generate images of galaxy pairs, using data from different
surveys.

"""

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

import attr

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import numpy as np

from retrying import retry

# =============================================================================
# METADATA
# =============================================================================

__version__ = "0.0.1"


# =============================================================================
# CONSTANTS
# =============================================================================

logger = logging.getLogger("scorpio")

#: Valid filters for the Sloan Digital Sky Survey (SDSS).
VALID_FILTERS_SDSS = ["u", "g", "r", "i", "z"]

#: Valid filters for the Two Micron All-Sky Survey (2MASS).
VALID_FILTERS_2MASS = ["-J", "-H", "-K"]

#: Valid filters for the Wide-field Infrared Survey Explorer (WISE).
VALID_FILTERS_WISE = [" 3.4", " 4.6", " 12", " 22"]

#: The default coordinate system.
DOWNLOAD_COORDINATE_SYSTEM = "J2000"

#: If a progress indicator will be showed when a file is downloaded.
DOWNLOAD_SHOW_PROGRESS = True

#: The default plot size for scorpio
PLOT_DEFAULT_SIZE = (8, 8)


# ============================================================================
# EXCEPTIONS
# ============================================================================


class NoFilterToStackError(ValueError):
    """Error generated when data of stack galaxies is empty."""

    pass


class NoSurveyInListToStackError(ValueError):
    """Error generated when survey of stack galaxies not in our list."""

    pass


# =============================================================================
# IMAGE RESPONSE CLASS
# =============================================================================


@attr.s(frozen=True)
class GPInteraction:
    """Representation of two interacting galaxies.

    Parameters
    ----------

    ra1 : float
        Right ascension of primary galaxy.
    dec1 : float
        Declination of primary galaxy.
    z1: float
        Redshift of secondary galaxy.
    ra2 : float
        Right ascension of secondary galaxy.
    dec2 : float
        Declination of secondary galaxy.
    z2 : float
        Redshift of secondary galaxy.
    resolution : int
        Size resolution value in pixels, by default it is 1000.
    survey : string
        Survey for query and download data, by default it is "SDSS".
    filters: list of strings
        Filters for respective survey, by default it is ["g", "i"] for "SDSS".

    Properties
    ----------
    mtx_ :
        The image of the two galaxies interacting.
    header_ :
        header information of image fits file.
    dist_physic_ :
        Estimation of physical distance to the observer in Mpc,
        from a given cosmology
    dist_pix_ :
        Physical distance as pixels in the image.
    pos1_ :
        RA, DEC, Z as pixels of the first galaxy.
    pos2_ :
        RA, DEC, Z as pixels of the second galaxy.


    """

    ra1 = attr.ib()
    dec1 = attr.ib()
    z1 = attr.ib()
    ra2 = attr.ib()
    dec2 = attr.ib()
    z2 = attr.ib()
    survey = attr.ib()
    resolution = attr.ib()
    cosmology = attr.ib()

    mtx_ = attr.ib(repr=False)
    header_ = attr.ib(repr=False)
    dist_physic_ = attr.ib(repr=False)
    dist_pix_ = attr.ib(repr=False)
    resolution_ = attr.ib(repr=False)
    pos1_ = attr.ib(repr=False)
    pos2_ = attr.ib(repr=False)

    def plot(
        self,
        ax=None,
        cmap="magma",
        center_color="cyan",
        cbar=True,
        center=True,
        scale=True,
        fmt=".4g",
        imshow_kws=None,
        cbar_kws=None,
        center_kws=None,
        center_text_kws=None,
        scalebar_kws=None,
        scalebar_text_kws=None,
        scalebar_bg_kws=None,
    ):
        r"""Receives data from other functions to generate and export a image.

        Parameters
        ----------
        ax : Axes
            Complete information of the properties to generate the image,
            by default it is None.
        cmap : ``str``, optional
            Name of the color map to be used
            (https://matplotlib.org/users/colormaps.html).
            If is None, the default colors of the matplotlib.pyplot.plot
            function is used, and if, and is a callable is used as
            colormap generator.
        center_color : ``str``, optional
            The color for the marker of the galaxy centers.
        cbar : ``bool``, optional.
            Whether to draw a colorbar.
        center : ``bool``, optional.
            Whether to draw the galaxy center markers.
        scale : ``bool``, optional.
            Whether to draw the scale bar.
        fmt : ``str``, optional
            String formatting code to use in the coordinates..
        imshow_kws : ``dict`` or ``None``, optional
             Keyword arguments for :meth:`matplotlib.pyplot.Axes.imshow`.
             This routine is the one that draws the image of the interaction.
        cbar_kws : ``dict`` or ``None``, optional
            Keyword arguments for :meth:`matplotlib.figure.Figure.colorbar`.
        center_kws : ``dict`` or ``None``, optional
            Keyword arguments for :meth:`matplotlib.pyplot.Axes.plot`.
            This routine is the one that draws the markers for the centers of the two galaxies.
        center_text_kws : ``dict`` or ``None``, optional
            Keyword arguments for :meth:`matplotlib.pyplot.Axes.text`.
            This routine is the one that draws the names of the markers of the galactic centers.
        scalebar_kws : ``dict`` or ``None``, optional
            Keyword arguments for :meth:`matplotlib.pyplot.Axes.hlines`.
            This routine is what draws the scale bar.
        scalebar_text_kws : ``dict`` or ``None``, optional
            Keyword arguments for :meth:`matplotlib.pyplot.Axes.text`.
            This routine is what draws the scale bar text.
        scalebar_bg_kws : ``dict`` or ``None``, optional
            Keyword arguments for :meth:`matplotlib.pyplot.Axes.broken_barh`.
            This routine is what draws the scale bar background color.

        Returns
        -------
        ``matplotlib.pyplot.Axis`` :
            The axis where the method draws.

        """
        if ax is None:
            ax, fig = plt.gca(), plt.gcf()
            fig.set_size_inches(PLOT_DEFAULT_SIZE)
        else:
            fig = ax.get_figure()

        # get al values from the instance for a more compact code
        survey = self.survey
        cosmology = self.cosmology
        ra1, dec1, z1 = [
            format(v, fmt) for v in (self.ra1, self.dec1, self.z1)
        ]
        ra2, dec2, z2 = [
            format(v, fmt) for v in (self.ra2, self.dec2, self.z2)
        ]

        final_image_a = np.copy(self.mtx_[0])
        plx = self.resolution_
        c1, c2 = self.pos1_, self.pos2_
        dis_c1_c2 = self.dist_pix_
        s_ab = self.dist_physic_

        xx = yy = plx / 2.0

        ax.axis([-xx * 0.8, xx * 0.8, -yy * 0.8, yy * 0.8])
        ax.xaxis.set_major_locator(ticker.NullLocator())
        ax.yaxis.set_major_locator(ticker.NullLocator())

        ax.set_title(
            f"Interaction - Survey: {survey} - " f"Cosmology: {cosmology.name}"
        )
        ax.set_xlabel("RA")
        ax.set_ylabel("Dec")

        # Normalization part 1
        max_value = np.max(final_image_a)
        min_value = np.min(final_image_a)

        limit = 0.0005 * max_value
        final_image_a[final_image_a <= limit] = limit

        # Normalization part 2
        max_value = np.max(final_image_a)
        min_value = np.min(final_image_a)

        # plot the image
        imshow_kws = {} if imshow_kws is None else dict(imshow_kws)
        imshow_kws["cmap"] = cmap
        imshow_kws.setdefault("extent", [-xx, xx, -yy, yy])
        imshow_kws.setdefault("norm", LogNorm(vmin=min_value, vmax=max_value))
        ax.imshow(final_image_a, **imshow_kws)

        # colorbar
        if cbar:
            cbar_kws = {} if cbar_kws is None else dict(cbar_kws)

            norm_color = mpl.colors.Normalize(vmin=min_value, vmax=max_value)
            scalar_mappable = mpl.cm.ScalarMappable(norm=norm_color, cmap=cmap)
            cax = inset_axes(
                ax,
                width=cbar_kws.pop("width", "5%"),
                height=cbar_kws.pop("height", "30%"),
                loc=cbar_kws.pop("loc", "lower right"),
            )

            cbar_kws.setdefault("mappable", scalar_mappable)
            cbar_kws.setdefault("cax", cax)
            cbar_kws.setdefault("orientation", "vertical")
            cb = fig.colorbar(ax=ax, **cbar_kws)
            cb.set_ticks([])

        # Scale bar background
        if scale:
            len_bar = 50.0 * dis_c1_c2 / s_ab

            scalebar_bg_kws = {} if scalebar_bg_kws is None else scalebar_bg_kws
            scalebar_bg_kws.setdefault(
                "xranges", [(-plx / 2.5, len_bar + plx / 7.5)]
            )
            scalebar_bg_kws.setdefault("yrange", (-plx / 2.5, plx / 7.5))
            scalebar_bg_kws.setdefault("facecolors", "w")
            ax.broken_barh(**scalebar_bg_kws)

            # Scale bar itself
            scalebar_kws = {} if scalebar_kws is None else dict(scalebar_kws)
            scalebar_kws.setdefault(
                "y",
                -plx / 2.72,
            )
            scalebar_kws.setdefault("xmin", -plx / 3.0)
            scalebar_kws.setdefault("xmax", -plx / 3.0 + len_bar)
            scalebar_kws.setdefault("color", "k")
            scalebar_kws.setdefault("linewidth", 3)
            ax.hlines(**scalebar_kws)

            # Scale bar text
            scalebar_text_kws = (
                {} if scalebar_text_kws is None else dict(scalebar_text_kws)
            )
            scalebar_text_kws.setdefault("x", -plx / 3.0)
            scalebar_text_kws.setdefault("y", -plx / 3.0)
            scalebar_text_kws.setdefault("s", "50 kpc")
            scalebar_text_kws.setdefault("fontsize", 20)
            ax.text(**scalebar_text_kws)

        # Markers of the galaxy centers
        if center:
            center_kws = {} if center_kws is None else dict(center_kws)
            center_kws["color"] = center_color
            center_kws.setdefault("marker", "o")
            center_kws.setdefault("markersize", 20)
            center_kws.setdefault("mew", 2)
            center_kws.setdefault("fillstyle", "none")

            label_g1 = f"G1: $RA={ra1}$, $Dec={dec1}$, $Z={z1}$"
            ax.plot(c1[0] - yy, -(c1[1] - xx), label=label_g1, **center_kws)

            label_g2 = f"G2: $RA={ra2}$, $Dec={dec2}$, $Z={z2}$"
            ax.plot(c2[0] - yy, -(c2[1] - xx), label=label_g2, **center_kws)

            center_text_kws = (
                {} if center_text_kws is None else dict(center_text_kws)
            )
            center_text_kws["color"] = center_color
            center_text_kws.setdefault("fontsize", 15)

            ax.text(c1[0] - yy + 20, -(c1[1] - xx), "G1", **center_text_kws)
            ax.text(c2[0] - yy + 20, -(c2[1] - xx), "G2", **center_text_kws)

            ax.legend(framealpha=1)

        return ax


# ============================================================================
# FUNCTIONS
# ============================================================================


@retry(stop_max_attempt_number=4)
def download_data(pos, survey, filters, plx, ff):
    """Download data fits from survey.

    Parameters
    ----------

    Returns
    -------

    """
    path = SkyView.get_images(
        position=pos,
        survey=survey + str(filters[ff]),
        radius=2 * apu.arcmin,
        pixels=(plx, plx),
        coordinates=DOWNLOAD_COORDINATE_SYSTEM,
        show_progress=DOWNLOAD_SHOW_PROGRESS,
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

    for default pixles: plx = 1000 and g, i filters.

    Parameters
    ----------
    ra1 : float
        Right ascension of primary galaxy.
    dec1 : float
        Declination of primary galaxy.
    z1: float
        Redshift of secondary galaxy.
    ra2 : float
        Right ascension of secondary galaxy.
    dec2 : float
        Declination of secondary galaxy.
    z2 : float
        Redshift of secondary galaxy.
    resolution : int
        Size resolution value in pixels, by default it is 1000.
    survey : string
        Survey for query and download data, by default it is "SDSS".
    filters: list of strings
        Filters for respective survey, by default it is ["g", "i"] for "SDSS".

    Returns
    -------
    g1g2 : array_like
        Array with the stacked information of the galaxies in the survey and
        the selected filters.
    header: string
        Header from the primary galaxy .fits file.
    plx : int
        Size resolution value in pixels.
    """

    plx = resolution

    glx1 = [ra1, dec1, z1]
    glx2 = [ra2, dec2, z2]

    glx_array = np.array([glx1, glx2])

    if survey not in ["SDSS", "2MASS", "WISE"]:
        raise NoSurveyInListToStackError(f"Survey '{survey}' not supported")

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

    Parameters
    ----------
    ra1 : float
        Right ascension of primary galaxy.
    dec1 : float
        Declination of primary galaxy.
    z1: float
        Redshift of secondary galaxy.
    ra2 : float
        Right ascension of secondary galaxy.
    dec2 : float
        Declination of secondary galaxy.
    z2 : float
        Redshift of secondary galaxy.
    header : string
        Header from the primary galaxy .fits file.
    survey : string
        Survey for query and download data, by default it is "SDSS".
    cosmology: astropy.cosmology.core.FlatLambdaCDM
        Instance of class ``astropy.cosmology.FLRW``,
        by default it is asc.Planck15.

    Returns
    -------
    s_ab : float
        physical distance between the pair in kpc.
    dis_c1_c2 : float
        physical distance to pixels in the image.
    c1, c2 : array_like
        Arrays with coordinates of both galaxies in pixels.
    """
    if cosmology not in [
        asc.WMAP5,
        asc.WMAP7,
        asc.WMAP9,
        asc.Planck13,
        asc.Planck15,
    ]:
        raise TypeError(f"cosmology `{cosmology}` not allowed")

    glx1 = [ra1, dec1, z1]
    glx2 = [ra2, dec2, z2]

    glx_array = np.array([glx1, glx2])
    z_glx = glx_array[:, 2]
    scale_factor = 1.0 / (1.0 + np.mean(z_glx))

    dist_comv = cosmology.comoving_distance(np.mean(z_glx)).value

    coord_a = SkyCoord(
        ra=glx_array[0, 0] * apu.deg, dec=glx_array[0, 1] * apu.deg
    )
    coord_b = SkyCoord(
        ra=glx_array[1, 0] * apu.deg, dec=glx_array[1, 1] * apu.deg
    )

    theta_rad = coord_a.separation(coord_b).rad
    s_ab = (dist_comv * theta_rad) * 1000.0 * scale_factor

    data_wcs = wcs.WCS(header)

    c1 = data_wcs.wcs_world2pix(glx_array[0, 0], glx_array[0, 1], 0)
    c2 = data_wcs.wcs_world2pix(glx_array[1, 0], glx_array[1, 1], 0)
    dis_c1_c2 = np.sqrt((c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2)

    return s_ab, dis_c1_c2, c1, c2


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

    Parameters
    ----------
    ra1 : float
        Right ascension of primary galaxy.
    dec1 : float
        Declination of primary galaxy.
    z1: float
        Redshift of secondary galaxy.
    ra2 : float
        Right ascension of secondary galaxy.
    dec2 : float
        Declination of secondary galaxy.
    z2 : float
        Redshift of secondary galaxy.
    survey : string
        Survey for query and download data, by default it is "SDSS".
    resolution : int
        Size resolution value in pixels, by default it is 1000.
    cosmology: astropy.cosmology.core.FlatLambdaCDM
        Instance of class astropy.cosmology.FLRW,
        by default it is asc.Planck15.

    Returns
    -------
    Image :
        An inmage of two interacting galaxies.

    """

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

    dist_physic, dist_pix, pos1, pos2 = distances(
        ra1, dec1, ra2, dec2, z1, z2, header, cosmology
    )

    gpi = GPInteraction(
        ra1=ra1,
        dec1=dec1,
        ra2=ra2,
        dec2=dec2,
        z1=z1,
        z2=z2,
        survey=survey,
        resolution=resolution,
        cosmology=cosmology,
        # properties
        mtx_=g1g2,
        header_=header,
        resolution_=plx,
        dist_physic_=dist_physic,
        dist_pix_=dist_pix,
        pos1_=pos1,
        pos2_=pos2,
    )

    return gpi
