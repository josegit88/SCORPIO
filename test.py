# -*- coding: utf-8 -*-
# Jose Benavides (jose.astroph@gmail.com)

# This file is part of the
#   SCORPIO Project (https://github.com/josegit88/SCORPIO).
# Copyright (c) 2020, Jose Benavides
# License: MIT
#   Full Text: https://github.com/josegit88/SCORPIO/blob/add-license-2/LICENSE

"""
Test about the scorpio functions
"""

import urllib

import astropy.cosmology as asc
from astropy.io import fits

from astroquery.skyview import SkyView

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm
from matplotlib.testing.decorators import check_figures_equal

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import numpy as np

import pytest

import scorpio


def test_download_and_stack(monkeypatch):
    [RA1, DEC1, Z1, RA2, DEC2, Z2] = [
        126.39162693999999,
        47.296980665521900,
        0.12573827000000001,
        126.38991429000001,
        47.305200665521902,
        0.12554201000000001,
    ]

    def mock_func_g1g2(position, survey, **kwargs):
        pos = [position.ra.value, position.dec.value]

        if pos == [RA1, DEC1] and survey == "SDSSg":
            path = [fits.open("test_data/SDSS_image_0_filter_g.fits")]
        elif pos == [RA1, DEC1] and survey == "SDSSi":
            path = [fits.open("test_data/SDSS_image_0_filter_i.fits")]
        elif pos == [RA2, DEC2] and survey == "SDSSg":
            path = [fits.open("test_data/SDSS_image_1_filter_g.fits")]
        elif pos == [RA2, DEC2] and survey == "SDSSi":
            path = [fits.open("test_data/SDSS_image_1_filter_i.fits")]

        return path

    monkeypatch.setattr(SkyView, "get_images", mock_func_g1g2)

    expected_g1 = np.array(
        [
            [-2.45361328e-02, -6.26220703e-02, -1.12152100e-02],
            [1.62239075e-02, 2.85937500e00, 1.81388855e-02],
            [-3.85131836e-02, 8.28247070e-02, 8.54492188e-04],
        ]
    )

    expected_g2 = np.array(
        [
            [1.13616943e-01, -1.66931152e-02, 4.88281250e-04],
            [1.61819458e-02, 2.73291016e00, 9.80529785e-02],
            [2.47344971e-02, -2.17952728e-02, 4.67758179e-02],
        ]
    )

    data_stack = scorpio.stack_pair(
        ra1=RA1, dec1=DEC1, ra2=RA2, dec2=DEC2, z1=Z1, z2=Z2, resolution=3
    )

    print(data_stack[0][0], data_stack[0][1])

    stack_g1 = data_stack[0][0]
    stack_g2 = data_stack[0][1]

    np.testing.assert_allclose(stack_g1, expected_g1, rtol=1e300)
    np.testing.assert_allclose(stack_g2, expected_g2, rtol=1e300)


# test for SDSS filters:
def test_download_invalid_filter_SDSS():
    [RA1, DEC1, Z1, RA2, DEC2, Z2] = [
        126.39162693999999,
        47.296980665521900,
        0.12573827000000001,
        126.38991429000001,
        47.305200665521902,
        0.12554201000000001,
    ]

    with pytest.raises(ValueError):
        scorpio.stack_pair(
            ra1=RA1,
            dec1=DEC1,
            ra2=RA2,
            dec2=DEC2,
            z1=Z1,
            z2=Z2,
            resolution=200,
            survey="SDSS",
            filters=["u", "t"],
        )


# test for 2MASS filters:
def test_download_invalid_filter_2MASS():
    [RA1, DEC1, Z1, RA2, DEC2, Z2] = [
        126.39162693999999,
        47.296980665521900,
        0.12573827000000001,
        126.38991429000001,
        47.305200665521902,
        0.12554201000000001,
    ]

    with pytest.raises(ValueError):
        scorpio.stack_pair(
            ra1=RA1,
            dec1=DEC1,
            ra2=RA2,
            dec2=DEC2,
            z1=Z1,
            z2=Z2,
            resolution=200,
            survey="2MASS",
            filters=["u", "t"],
        )


# test for WISE filters:
def test_download_invalid_filter_WISE():
    [RA1, DEC1, Z1, RA2, DEC2, Z2] = [
        126.39162693999999,
        47.296980665521900,
        0.12573827000000001,
        126.38991429000001,
        47.305200665521902,
        0.12554201000000001,
    ]

    with pytest.raises(ValueError):
        scorpio.stack_pair(
            ra1=RA1,
            dec1=DEC1,
            ra2=RA2,
            dec2=DEC2,
            z1=Z1,
            z2=Z2,
            resolution=200,
            survey="WISE",
            filters=["u", "t"],
        )


# new:
def test_download_invalid_survey():
    [RA1, DEC1, Z1, RA2, DEC2, Z2] = [
        126.39162693999999,
        47.296980665521900,
        0.12573827000000001,
        126.38991429000001,
        47.305200665521902,
        0.12554201000000001,
    ]

    with pytest.raises(ValueError):
        scorpio.stack_pair(
            ra1=RA1,
            dec1=DEC1,
            ra2=RA2,
            dec2=DEC2,
            z1=Z1,
            z2=Z2,
            resolution=200,
            survey="JABB",
        )


def test_stack_code_error(monkeypatch):
    def mock_func(*args, **kwargs):
        raise urllib.error.HTTPError(
            url="http://from.mock",
            code=500,
            msg="from mock",
            hdrs=None,
            fp=None,
        )

    monkeypatch.setattr(SkyView, "get_images", mock_func)

    [RA1, DEC1, Z1, RA2, DEC2, Z2] = [
        229.38791793999997,
        -15.1525941155219059,
        0.12589231000000001,
        229.38705890000003,
        -15.1513408255219060,
        0.12634666999999999,
    ]

    with pytest.raises(urllib.error.HTTPError):
        scorpio.stack_pair(
            ra1=RA1,
            dec1=DEC1,
            ra2=RA2,
            dec2=DEC2,
            z1=Z1,
            z2=Z2,
            resolution=450,
            filters=["g"],
        )


def test_distances_error_Cosmology():
    data_imagen = fits.open("test_data/SDSS_image_0_filter_g.fits")
    header = data_imagen[0].header
    COSMOLOGY = {"H0": 70, "Om0": 0.3, "Ode0": 0.7}

    [RA1, DEC1, Z1, RA2, DEC2, Z2] = [
        126.39162693999999,
        47.296980665521900,
        0.12573827000000001,
        126.38991429000001,
        47.305200665521902,
        0.12554201000000001,
    ]

    with pytest.raises(TypeError):
        scorpio.distances(
            ra1=RA1,
            dec1=DEC1,
            ra2=RA2,
            dec2=DEC2,
            z1=Z1,
            z2=Z2,
            header=header,
            cosmology=COSMOLOGY,
        )


def test_distances_physical():
    data_imagen = fits.open("test_data/SDSS_image_0_filter_g.fits")
    header = data_imagen[0].header
    COSMOLOGY = asc.Planck15

    [RA1, DEC1, Z1, RA2, DEC2, Z2] = [
        126.39162693999999,
        47.296980665521900,
        0.12573827000000001,
        126.38991429000001,
        47.305200665521902,
        0.12554201000000001,
    ]

    data_distances = scorpio.distances(
        ra1=RA1,
        dec1=DEC1,
        ra2=RA2,
        dec2=DEC2,
        z1=Z1,
        z2=Z2,
        header=header,
        cosmology=COSMOLOGY,
    )

    physical_dist = data_distances[0]
    expected_dist = 78.169
    np.testing.assert_allclose(physical_dist, expected_dist, rtol=1e-2)


def test_distances_pixels():
    data_imagen = fits.open("test_data/SDSS_image_0_filter_g.fits")
    header = data_imagen[0].header
    COSMOLOGY = asc.Planck15

    [RA1, DEC1, Z1, RA2, DEC2, Z2] = [
        126.39162693999999,
        47.296980665521900,
        0.12573827000000001,
        126.38991429000001,
        47.305200665521902,
        0.12554201000000001,
    ]

    data_distances = scorpio.distances(
        ra1=RA1,
        dec1=DEC1,
        ra2=RA2,
        dec2=DEC2,
        z1=Z1,
        z2=Z2,
        header=header,
        cosmology=COSMOLOGY,
    )

    pixel_dist = data_distances[1]
    expected_dist = 0.747148061
    np.testing.assert_allclose(pixel_dist, expected_dist, rtol=1e-5)


# test kwargs:
def test_plot_kwargs():
    [RA1, DEC1, Z1, RA2, DEC2, Z2] = [
        126.39162693999999,
        47.296980665521900,
        0.12573827000000001,
        126.38991429000001,
        47.305200665521902,
        0.12554201000000001,
    ]

    test_img = scorpio.gpair(
        ra1=RA1,
        dec1=DEC1,
        ra2=RA2,
        dec2=DEC2,
        z1=Z1,
        z2=Z2,
        survey="2MASS",
        resolution=500,
    )

    with pytest.raises(AttributeError):
        test_img.plot("./dir_test_images")

    with pytest.raises(AttributeError):
        test_img.plot("y")

    with pytest.raises(AttributeError):
        test_img.plot("img_test.png")


def test_plot_ax():
    [RA1, DEC1, Z1, RA2, DEC2, Z2] = [
        126.39162693999999,
        47.296980665521900,
        0.12573827000000001,
        126.38991429000001,
        47.305200665521902,
        0.12554201000000001,
    ]

    test_img = scorpio.gpair(
        ra1=RA1,
        dec1=DEC1,
        ra2=RA2,
        dec2=DEC2,
        z1=Z1,
        z2=Z2,
        survey="2MASS",
        resolution=500,
    )

    with pytest.raises(AttributeError):
        test_img.plot(
            ax=plt.subplots(figsize=(4, 8)),
        )


@check_figures_equal(extensions=["png"])
def test_equal_download_and_generate_plots(fig_test, fig_ref):
    [RA1, DEC1, Z1, RA2, DEC2, Z2] = [
        126.39162693999999,
        47.296980665521900,
        0.12573827000000001,
        126.38991429000001,
        47.305200665521902,
        0.12554201000000001,
    ]

    data_img1 = scorpio.gpair(
        ra1=RA1, dec1=DEC1, ra2=RA2, dec2=DEC2, z1=Z1, z2=Z2
    )
    data_img2 = scorpio.gpair(
        ra1=RA1,
        dec1=DEC1,
        ra2=RA2,
        dec2=DEC2,
        z1=Z1,
        z2=Z2,
    )

    # test plot 1:
    test_ax1 = fig_test.subplots()
    ax1 = data_img1.plot(ax=test_ax1)
    ax1.set_ylabel("DEC")
    ax1.set_xlabel("RA")

    # test plot 2
    test_ax2 = fig_ref.subplots()
    ax2 = data_img2.plot(ax=test_ax2)
    ax2.set_ylabel("DEC")
    ax2.set_xlabel("RA")


@check_figures_equal(extensions=["png"])
def test_equal_plots(fig_test, fig_ref):
    [RA1, DEC1, Z1, RA2, DEC2, Z2] = [
        126.39162693999999,
        47.296980665521900,
        0.12573827000000001,
        126.38991429000001,
        47.305200665521902,
        0.12554201000000001,
    ]

    data_img = scorpio.gpair(
        ra1=RA1,
        dec1=DEC1,
        ra2=RA2,
        dec2=DEC2,
        z1=Z1,
        z2=Z2,
        survey="2MASS",
        resolution=500,
    )

    # fig test
    test_ax = fig_test.subplots()
    ax1 = data_img.plot(ax=test_ax)
    ax1.set_ylabel("DEC")
    ax1.set_xlabel("RA")

    # fig expect
    final_imageA = data_img.matriz
    final_imageA = final_imageA[0]
    plx = data_img.resolution
    c1 = data_img.pos1
    c2 = data_img.pos2
    dis_c1_c2 = data_img.dist_pix
    s_AB = data_img.dist_physic

    expect_ax = fig_ref.subplots()
    fig = plt.gcf()

    xx = plx / 2.0
    yy = plx / 2.0

    expect_ax.axis([-xx * 0.8, xx * 0.8, -yy * 0.8, yy * 0.8])
    expect_ax.xaxis.set_major_locator(ticker.NullLocator())
    expect_ax.yaxis.set_major_locator(ticker.NullLocator())

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
    Cmap = mpl.cm.ScalarMappable(norm=Norm, cmap=mpl.cm.inferno)
    axins1 = inset_axes(
        expect_ax,
        width="5%",
        height="30%",
        loc="lower right",
    )
    cb = fig.colorbar(Cmap, ax=expect_ax, cax=axins1, orientation="vertical")
    cb.set_ticks([])
    expect_ax.imshow(
        final_imageA,
        extent=extent,
        cmap="inferno",
        norm=LogNorm(vmin=min_value, vmax=max_value),
    )
    expect_ax.plot(
        c1[0] - yy,
        -(c1[1] - xx),
        color="cyan",
        marker="o",
        markersize=20,
        mew=2,
        fillstyle="none",
    )
    expect_ax.plot(
        c2[0] - yy,
        -(c2[1] - xx),
        color="cyan",
        marker="o",
        markersize=20,
        mew=2,
        fillstyle="none",
    )

    len_bar = 50.0 * dis_c1_c2 / s_AB
    expect_ax.broken_barh(
        [(-plx / 2.5, len_bar + plx / 7.5)],
        (-plx / 2.5, plx / 7.5),
        facecolors="w",
    )
    expect_ax.hlines(
        y=-plx / 2.72,
        xmin=-plx / 3.0,
        xmax=-plx / 3.0 + len_bar,
        color="k",
        linewidth=3,
    )
    expect_ax.text(-plx / 3.0, -plx / 3.0, "50 kpc", fontsize=20, color="k")
    expect_ax.set_ylabel("DEC")
    expect_ax.set_xlabel("RA")
    # ----------------------------------
