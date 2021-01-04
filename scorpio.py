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
import urllib

import astropy.cosmology as asc
from astropy import units as apu
from astropy import wcs
from astropy.coordinates import SkyCoord

from astroquery.skyview import SkyView

import attr

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import numpy as np

from retrying import retry

import seaborn as sns


# =============================================================================
# METADATA
# =============================================================================

__version__ = "0.2"


# =============================================================================
# CONSTANTS
# =============================================================================

logger = logging.getLogger("scorpio")

#: SURVEYS
SURVEYS = {
    "SDSS": {"default": ["g", "i"], "filters": ["u", "g", "r", "i", "z"]},
    "2MASS": {"default": ["-H", "-K"], "filters": ["-J", "-H", "-K"]},
    "WISE": {
        "default": [" 4.6", " 12"],
        "filters": [" 3.4", " 4.6", " 12", " 22"],
    },
}

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

    Attributes
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
    pos1_ = attr.ib(repr=False)
    pos2_ = attr.ib(repr=False)

    def plot(
        self,
        ax=None,
        cmap="magma",
        center_color="cyan",
        center=True,
        scale=True,
        fmt=".4g",
        llimit=None,
        heatmap_kws=None,
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
        center : ``bool``, optional.
            Whether to draw the galaxy center markers.
        scale : ``bool``, optional.
            Whether to draw the scale bar.
        fmt : ``str``, optional
            String formatting code to use in the coordinates.
        llimit : ``float``, optional.
            The lowest value to show in the plot. If some value is <= llimit
            then it was remplaced with llimit. By default this value is setted
            as 0.0005 times the maximun value of the image.
        heatmap_kws : ``dict`` or ``None``, optional
             Keyword arguments for :meth:`seaborn.heatmap`.
             This routine is the one that draws the image of the interaction.
        cbar_kws : ``dict`` or ``None``, optional
            Keyword arguments for :meth:`matplotlib.figure.Figure.colorbar`.
        center_kws : ``dict`` or ``None``, optional
            Keyword arguments for :meth:`matplotlib.pyplot.Axes.plot`.
            This routine is the one that draws the markers for the centers of
            the two galaxies.
        center_text_kws : ``dict`` or ``None``, optional
            Keyword arguments for :meth:`matplotlib.pyplot.Axes.text`.
            This routine is the one that draws the names of the markers of the
            galactic centers.
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

        # get al values from the instance for a more compact code
        survey = self.survey
        cosmology = self.cosmology

        final_image_a = np.copy(self.mtx_[0])
        plx = self.resolution
        c1, c2 = self.pos1_, self.pos2_
        dis_c1_c2 = self.dist_pix_
        s_ab = self.dist_physic_.value

        # Normalization part 1
        max_value = np.max(final_image_a)
        min_value = np.min(final_image_a)

        # Normalization part 2 (remove negatives)
        llimit = (0.0005 * max_value) if llimit is None else llimit
        final_image_a[final_image_a <= llimit] = llimit

        max_value = np.max(final_image_a)
        min_value = np.min(final_image_a)

        # plot the image
        heatmap_kws = {} if heatmap_kws is None else dict(heatmap_kws)
        heatmap_kws["cmap"] = cmap
        heatmap_kws.setdefault("square", True)
        heatmap_kws.setdefault("cbar", False)
        heatmap_kws.setdefault("norm", LogNorm(vmin=min_value, vmax=max_value))
        heatmap_kws.setdefault("xticklabels", False)
        heatmap_kws.setdefault("yticklabels", False)
        sns.heatmap(final_image_a, ax=ax, **heatmap_kws)

        # Markers of the galaxy centers
        if center:
            ra1, dec1, z1 = [
                format(v, fmt) for v in (self.ra1, self.dec1, self.z1)
            ]
            ra2, dec2, z2 = [
                format(v, fmt) for v in (self.ra2, self.dec2, self.z2)
            ]

            center_kws = {} if center_kws is None else dict(center_kws)
            center_kws["color"] = center_color
            center_kws.setdefault("marker", "o")
            center_kws.setdefault("markersize", 20)
            center_kws.setdefault("mew", 2)
            center_kws.setdefault("fillstyle", "none")

            label_g1 = f"G1: $RA={ra1}$, $Dec={dec1}$, $Z={z1}$"
            ax.plot(c1[0], c1[1], label=label_g1, **center_kws)

            label_g2 = f"G2: $RA={ra2}$, $Dec={dec2}$, $Z={z2}$"
            ax.plot(c2[0], c2[1], label=label_g2, **center_kws)

            center_text_kws = (
                {} if center_text_kws is None else dict(center_text_kws)
            )
            center_text_kws["color"] = center_color
            center_text_kws.setdefault("fontsize", 15)

            offset = center_kws.get("markersize") or 0
            ax.text(c1[0] + offset, c1[1], "G1", **center_text_kws)
            ax.text(c2[0] + offset, c2[1], "G2", **center_text_kws)

            ax.legend(framealpha=1)

        if scale:
            len_bar = 50.0 * dis_c1_c2 / s_ab

            # Scale bar kwargs text
            scalebar_text_kws = (
                {} if scalebar_text_kws is None else dict(scalebar_text_kws)
            )
            scale_fontsize = scalebar_text_kws.setdefault("fontsize", 15)
            scalebar_text_kws.setdefault("x", scale_fontsize)
            scalebar_text_kws.setdefault("y", plx - scale_fontsize * 2)
            scalebar_text_kws.setdefault("s", "50 kpc")
            scalebar_text_kws.setdefault("color", "k")

            # Scale bar kwargs itself
            scalebar_kws = {} if scalebar_kws is None else dict(scalebar_kws)
            scalebar_kws.setdefault("y", plx - scale_fontsize * 1.1)
            scalebar_kws.setdefault("xmin", scale_fontsize)
            scalebar_kws.setdefault("xmax", len_bar)
            scalebar_kws.setdefault("color", "k")
            scalebar_kws.setdefault("linewidth", 3)

            # Scale bar background kwargs
            scalebar_bg_kws = (
                {} if scalebar_bg_kws is None else scalebar_bg_kws
            )
            scalebar_bg_kws.setdefault("xranges", [(0, len_bar * 1.2)])
            scalebar_bg_kws.setdefault(
                "yrange", (plx * 0.9 - scale_fontsize, plx)
            )
            scalebar_bg_kws.setdefault("facecolors", "w")

            ax.broken_barh(**scalebar_bg_kws)
            ax.hlines(**scalebar_kws)
            ax.text(**scalebar_text_kws)

        ax.set_title(
            f"Interaction - Survey: {survey} - " f"Cosmology: {cosmology.name}"
        )
        ax.set_xlabel("RA")
        ax.set_ylabel("Dec")

        ax.patch.set_edgecolor("black")
        ax.patch.set_linewidth("3")

        return ax


# ============================================================================
# FUNCTIONS
# ============================================================================


@retry(stop_max_attempt_number=4)
def download_data(pos, survey, filters, plx):
    """Proxy to astroquery SkyviewService.

    This functions is mostly for internal use.

    Parameters
    ----------
    pos: SkyCoord
        Determines the center of the field to be retrieved.
    survey: str
        The data to download the survey.
    filters: list
        Specific astronomical filters of the data.
    plx:
        Selects the pixel dimensions of the image to be produced.

    Returns
    -------
    A list of `~astropy.io.fits.HDUList` objects.

    """
    path = SkyView.get_images(
        position=pos,
        survey=survey + str(filters),
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
    # survey validation
    survey_conf = SURVEYS.get(survey)
    if survey_conf is None:
        raise NoSurveyInListToStackError(f"Survey '{survey}' not supported")

    filters = survey_conf["default"] if filters is None else filters
    for ff in filters:
        if ff not in survey_conf["filters"]:
            raise NoFilterToStackError(f"{ff} is not a valid filter")

    # preprocess
    plx = resolution

    glx1 = [ra1, dec1, z1]
    glx2 = [ra2, dec2, z2]

    glx_array = np.array([glx1, glx2])

    # ------- download images data fits in diferent filters: --------
    stamp0 = None
    g1g2 = [np.zeros(shape=(plx, plx)), np.zeros(shape=(plx, plx))]

    for ff in range(len(filters)):
        for ii in range(len(g1g2)):
            pos = SkyCoord(
                ra=glx_array[ii, 0] * apu.degree,
                dec=glx_array[ii, 1] * apu.degree,
            )
            try:
                stamp = download_data(
                    pos=pos, survey=survey, filters=filters[ff], plx=plx
                )
                if ii == 0:
                    stamp0 = stamp
            except urllib.error.HTTPError as err:
                if err.code != 404:
                    raise err
                logger.warning(f"No data filter '{filters[ff]}' in gal {ii}")
            else:
                g1g2[ii] += stamp[0][0].data

    if np.all(g1g2[0] == 0):
        raise NoFilterToStackError("Empty array for galaxy1")

    header = stamp0[0][0].header
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
    survey : string, optional (default='SDSS')
        Survey for query and download data.
    resolution : int, optional (default=1000)
        Size resolution value in pixels.
    cosmology: astropy.cosmology.core.FlatLambdaCDM, optional
        Instance of class ``astropy.cosmology.FLRW``.
        default `asc.Planck15`.

    Returns
    -------
    GPInteraction :
        An object with information about the two interacting galaxies.

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
        dist_physic_=dist_physic * apu.Mpc,
        dist_pix_=dist_pix,
        pos1_=pos1,
        pos2_=pos2,
    )

    return gpi
