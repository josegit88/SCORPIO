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

import os
import pathlib
import urllib

import astropy.cosmology as apc
from astropy.io import fits

from astroquery.skyview import SkyView

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.testing.decorators import check_figures_equal

import numpy as np

import pytest

import scorpio

import seaborn as sns

# =============================================================================
# CONSTANTS
# =============================================================================

PATH = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))

TEST_DATA = PATH / "test_data"

# =============================================================================
# FIXTURE
# =============================================================================


@pytest.fixture(scope="session")
def mget_images():
    def _mock(position, survey, **kwargs):
        pos = str([position.ra.value, position.dec.value])

        call_str = (
            f"SkyView.get_images(position={pos}, survey='{survey}', "
            f"radius='{kwargs['radius']}', pixels={kwargs['pixels']}, "
            f"coordinates='{kwargs['coordinates']}', "
            f"show_progress={kwargs['show_progress']})"
        )
        full_path = TEST_DATA / "mget_images" / call_str / "return.fits"

        return [fits.open(full_path)]

    return _mock


# =============================================================================
# TESTS
# =============================================================================


def test_download_and_stack_data(mget_images, monkeypatch):

    ra1, dec1, ra2, dec2 = [
        126.39162693999999,
        47.296980665521900,
        126.38991429000001,
        47.305200665521902,
    ]

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

    monkeypatch.setattr(SkyView, "get_images", mget_images)

    data_stack = scorpio.stack_pair(
        ra1=ra1, dec1=dec1, ra2=ra2, dec2=dec2, resolution=3
    )

    stack_g1 = data_stack[0][0]
    stack_g2 = data_stack[0][1]

    np.testing.assert_allclose(stack_g1, expected_g1)
    np.testing.assert_allclose(stack_g2, expected_g2)


# test for SDSS filters:
def test_download_invalid_filter_SDSS():
    [ra1, dec1, ra2, dec2] = [
        126.39162693999999,
        47.296980665521900,
        126.38991429000001,
        47.305200665521902,
    ]

    with pytest.raises(ValueError):
        scorpio.stack_pair(
            ra1=ra1,
            dec1=dec1,
            ra2=ra2,
            dec2=dec2,
            resolution=200,
            survey="SDSS",
            filters=["u", "t"],
        )


# test for 2MASS filters:
def test_download_invalid_filter_2MASS():
    [ra1, dec1, ra2, dec2] = [
        126.39162693999999,
        47.296980665521900,
        126.38991429000001,
        47.305200665521902,
    ]

    with pytest.raises(ValueError):
        scorpio.stack_pair(
            ra1=ra1,
            dec1=dec1,
            ra2=ra2,
            dec2=dec2,
            resolution=200,
            survey="2MASS",
            filters=["u", "t"],
        )


# test for WISE filters:
def test_download_invalid_filter_WISE():
    [ra1, dec1, ra2, dec2] = [
        126.39162693999999,
        47.296980665521900,
        126.38991429000001,
        47.305200665521902,
    ]

    with pytest.raises(ValueError):
        scorpio.stack_pair(
            ra1=ra1,
            dec1=dec1,
            ra2=ra2,
            dec2=dec2,
            resolution=200,
            survey="WISE",
            filters=["u", "t"],
        )


def test_download_invalid_survey():
    [ra1, dec1, ra2, dec2] = [
        126.39162693999999,
        47.296980665521900,
        126.38991429000001,
        47.305200665521902,
    ]

    with pytest.raises(ValueError):
        scorpio.stack_pair(
            ra1=ra1,
            dec1=dec1,
            ra2=ra2,
            dec2=dec2,
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

    [ra1, dec1, ra2, dec2] = [
        229.38791793999997,
        -15.1525941155219059,
        229.38705890000003,
        -15.1513408255219060,
    ]

    with pytest.raises(urllib.error.HTTPError):
        scorpio.stack_pair(
            ra1=ra1,
            dec1=dec1,
            ra2=ra2,
            dec2=dec2,
            resolution=450,
            filters=["g"],
        )


def test_distances_error_Cosmology():
    image = fits.open(TEST_DATA / "SDSS_g.fits")
    header = image[0].header
    cosmology = {"H0": 70, "Om0": 0.3, "Ode0": 0.7}

    [ra1, dec1, z1, ra2, dec2, z2] = [
        126.39162693999999,
        47.296980665521900,
        0.12573827000000001,
        126.38991429000001,
        47.305200665521902,
        0.12554201000000001,
    ]

    with pytest.raises(TypeError):
        scorpio.distances(
            ra1=ra1,
            dec1=dec1,
            ra2=ra2,
            dec2=dec2,
            z1=z1,
            z2=z2,
            header=header,
            cosmology=cosmology,
        )


def test_distances_physical_units():
    image = fits.open(TEST_DATA / "SDSS_g.fits")
    header = image[0].header
    cosmology = apc.Planck15

    [ra1, dec1, z1, ra2, dec2, z2] = [
        126.39162693999999,
        47.296980665521900,
        0.12573827000000001,
        126.38991429000001,
        47.305200665521902,
        0.12554201000000001,
    ]

    data_distances = scorpio.distances(
        ra1=ra1,
        dec1=dec1,
        ra2=ra2,
        dec2=dec2,
        z1=z1,
        z2=z2,
        header=header,
        cosmology=cosmology,
    )

    physical_dist = data_distances[0]
    expected_dist = 69.444
    np.testing.assert_allclose(physical_dist, expected_dist, rtol=1e-2)


def test_distances_in_pixels():
    image = fits.open(TEST_DATA / "SDSS_g.fits")
    header = image[0].header
    cosmology = apc.Planck15

    [ra1, dec1, z1, ra2, dec2, z2] = [
        126.39162693999999,
        47.296980665521900,
        0.12573827000000001,
        126.38991429000001,
        47.305200665521902,
        0.12554201000000001,
    ]

    data_distances = scorpio.distances(
        ra1=ra1,
        dec1=dec1,
        ra2=ra2,
        dec2=dec2,
        z1=z1,
        z2=z2,
        header=header,
        cosmology=cosmology,
    )

    pixel_dist = data_distances[1]
    expected_dist = 0.747148061
    np.testing.assert_allclose(pixel_dist, expected_dist, rtol=1e-5)


def test_plot_size_fig_axes(mget_images, monkeypatch):

    monkeypatch.setattr(SkyView, "get_images", mget_images)

    ra1, dec1, z1, ra2, dec2, z2 = [
        126.39162693999999,
        47.296980665521900,
        0.12573827000000001,
        126.38991429000001,
        47.305200665521902,
        0.12554201000000001,
    ]

    test_img = scorpio.gpair(
        ra1=ra1,
        dec1=dec1,
        ra2=ra2,
        dec2=dec2,
        z1=z1,
        z2=z2,
        resolution=3,
    )

    with pytest.raises(AttributeError):
        test_img.plot(
            ax=plt.subplots(figsize=(4, 8)),
        )


@check_figures_equal(extensions=["png"])
def test_download_and_generate_equal_plots(
    fig_test, fig_ref, mget_images, monkeypatch
):

    monkeypatch.setattr(SkyView, "get_images", mget_images)

    ra1, dec1, z1, ra2, dec2, z2 = [
        126.39162693999999,
        47.296980665521900,
        0.12573827000000001,
        126.38991429000001,
        47.305200665521902,
        0.12554201000000001,
    ]

    data_img1 = scorpio.gpair(
        ra1=ra1,
        dec1=dec1,
        ra2=ra2,
        dec2=dec2,
        z1=z1,
        z2=z2,
        resolution=3,
    )

    data_img2 = scorpio.gpair(
        ra1=ra1,
        dec1=dec1,
        ra2=ra2,
        dec2=dec2,
        z1=z1,
        z2=z2,
        resolution=3,
    )

    # test plot 1:
    test_ax1 = fig_test.subplots()
    data_img1.plot(ax=test_ax1)

    # test plot 2
    test_ax2 = fig_ref.subplots()
    data_img2.plot(ax=test_ax2)


@check_figures_equal()
def test_compare_plots_generation_methods(
    fig_test, fig_ref, mget_images, monkeypatch
):

    monkeypatch.setattr(SkyView, "get_images", mget_images)

    ra1, dec1, z1, ra2, dec2, z2 = [
        126.39162693999999,
        47.296980665521900,
        0.12573827000000001,
        126.38991429000001,
        47.305200665521902,
        0.12554201000000001,
    ]

    data_img = scorpio.gpair(
        ra1=ra1,
        dec1=dec1,
        ra2=ra2,
        dec2=dec2,
        z1=z1,
        z2=z2,
        resolution=3,
    )

    # >>>>>>>>> fig test
    test_ax = fig_test.subplots()
    data_img.plot(ax=test_ax)

    # >>>>>>>>> fig expect
    expect_ax = fig_ref.subplots()
    survey = data_img.survey
    cosmology = data_img.cosmology

    final_image_a = np.copy(data_img.mtx_[0])
    plx = data_img.resolution
    c1, c2 = data_img.pos1_, data_img.pos2_
    dis_c1_c2 = data_img.dist_pix_
    s_ab = data_img.dist_physic_.value

    # Normalization part 1
    max_value = np.max(final_image_a)
    min_value = np.min(final_image_a)

    # Normalization part 2 (remove negatives)
    llimit = 0.0005 * max_value
    final_image_a[final_image_a <= llimit] = llimit

    max_value = np.max(final_image_a)
    min_value = np.min(final_image_a)

    # plot the image
    sns.heatmap(
        final_image_a,
        ax=expect_ax,
        cmap="magma",
        square=True,
        cbar=False,
        norm=LogNorm(vmin=min_value, vmax=max_value),
        xticklabels=False,
        yticklabels=False,
    )

    # Markers of the galaxy centers
    fmt = ".4g"
    ra1, dec1, z1 = [
        format(v, fmt) for v in (data_img.ra1, data_img.dec1, data_img.z1)
    ]
    ra2, dec2, z2 = [
        format(v, fmt) for v in (data_img.ra2, data_img.dec2, data_img.z2)
    ]

    center_kws = {
        "color": "cyan",
        "marker": "o",
        "markersize": 20,
        "mew": 2,
        "fillstyle": "none",
    }
    label_g1 = f"G1: $RA={ra1}$, $Dec={dec1}$, $Z={z1}$"
    expect_ax.plot(c1[0], c1[1], label=label_g1, **center_kws)

    label_g2 = f"G2: $RA={ra2}$, $Dec={dec2}$, $Z={z2}$"
    expect_ax.plot(c2[0], c2[1], label=label_g2, **center_kws)

    center_text_kws = {"color": "cyan", "fontsize": 15}

    offset = 20
    expect_ax.text(c1[0] + offset, c1[1], "G1", **center_text_kws)
    expect_ax.text(c2[0] + offset, c2[1], "G2", **center_text_kws)
    expect_ax.legend(framealpha=1)

    # scalebar
    len_bar = 50.0 * dis_c1_c2 / s_ab

    # Scale bar background
    expect_ax.broken_barh(  # expect_ax.broken_barh(**scalebar_bg_kws)
        xranges=[(0, len_bar * 1.2)],
        yrange=(plx * 0.9 - 15, plx),
        facecolors="w",
    )

    # Scale bar kwargs itself
    expect_ax.hlines(
        y=plx - 15 * 1.1,
        xmin=15,
        xmax=len_bar,
        color="k",
        linewidth=3,
    )

    expect_ax.text(
        fontsize=15,
        x=15,
        y=plx - 15 * 2,
        s="50 kpc",
        color="k",
    )

    expect_ax.set_title(
        f"Interaction - Survey: {survey} - " f"Cosmology: {cosmology.name}"
    )
    expect_ax.set_xlabel("RA")
    expect_ax.set_ylabel("Dec")

    expect_ax.patch.set_edgecolor("black")
    expect_ax.patch.set_linewidth("3")
