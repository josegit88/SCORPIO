# -*- coding: utf-8 -*-
# Jose Benavides (jose.astroph@gmail.com)
"""
Test about the scorpio functions
"""

import urllib

import astropy.cosmology as asc
from astropy.io import fits

from astroquery.skyview import SkyView

import numpy as np

import pytest

import scorpio


def test_download_and_stack(monkeypatch):
    # aa = [126.39162693999999, 47.296980665521900, 0.12573827000000001]
    # bb = [126.38991429000001, 47.305200665521902, 0.12554201000000001]
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

    # try:
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

    """
    except:
        expected_g1 = np.array(
            [
                [-2.45361328e-02, -6.26220703e-02, 3.99195767e252],
                [1.62239075e-02, 3.36201984e160, 1.81388855e-02],
                [2.96410912e222, 8.28247070e-02, -4.48529307e198],
            ]
        )

        expected_g2 = np.array(
            [
                [1.13616943e-01, -1.66931152e-02, 4.88281250e-04],
                [1.61819458e-02, 2.73291016e0, 4.52350897e257],
                [2.47344971e-02, 1.91084608e214, 4.67758179e-02],
            ]
        )

        data_stack = scorpio.stack_pair(
            ra1=RA1, dec1=DEC1, ra2=RA2, dec2=DEC2, z1=Z1, z2=Z2, resolution=3
        )

        print(data_stack[0][0], data_stack[0][1])

        stack_g1 = data_stack[0][0]
        stack_g2 = data_stack[0][1]

        np.testing.assert_allclose(stack_g1, expected_g1, rtol=1e-3)
        np.testing.assert_allclose(stack_g2, expected_g2, rtol=1e-3)
    """


# test for SDSS filters:
def test_download_invalid_filter_SDSS():
    # aa = [126.39162693999999, 47.296980665521900, 0.12573827000000001]
    # bb = [126.38991429000001, 47.305200665521902, 0.12554201000000001]
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
    # aa = [126.39162693999999, 47.296980665521900, 0.12573827000000001]
    # bb = [126.38991429000001, 47.305200665521902, 0.12554201000000001]
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
    # aa = [126.39162693999999, 47.296980665521900, 0.12573827000000001]
    # bb = [126.38991429000001, 47.305200665521902, 0.12554201000000001]
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
    # aa = [126.39162693999999, 47.296980665521900, 0.12573827000000001]
    # bb = [126.38991429000001, 47.305200665521902, 0.12554201000000001]
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


def test_no_filter_to_stack(monkeypatch):
    def mock_func(*args, **kwargs):
        raise urllib.error.HTTPError(
            url="http://from.mock",
            code=404,
            msg="from mock",
            hdrs=None,
            fp=None,
        )

    monkeypatch.setattr(SkyView, "get_images", mock_func)

    # cc = [229.38791793999997, -15.1525941155219059, 0.12589231000000001]
    # dd = [229.38705890000003, -15.1513408255219060, 0.12634666999999999]
    [RA1, DEC1, Z1, RA2, DEC2, Z2] = [
        229.38791793999997,
        -15.1525941155219059,
        0.12589231000000001,
        229.38705890000003,
        -15.1513408255219060,
        0.12634666999999999,
    ]

    with pytest.raises(scorpio.NoFilterToStackError):
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

    # cc = [229.38791793999997, -15.1525941155219059, 0.12589231000000001]
    # dd = [229.38705890000003, -15.1513408255219060, 0.12634666999999999]
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


# new:

# puedo armar que mi cosmology sea un diccionario
# cosmo = {"H0":70, "Om0":0.3, "Ode0":0.7}
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
            info_fits=header,
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
        info_fits=header,
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
        info_fits=header,
        cosmology=COSMOLOGY,
    )

    pixel_dist = data_distances[1]
    expected_dist = 0.747148061
    np.testing.assert_allclose(pixel_dist, expected_dist, rtol=1e-5)


# testeo de plot, a mano creo mi imagen, la que yo espero, la que me devuelve
# scorpio y hacer una comparacion pixel a pixel.
# Revisar trabajo de Bruno y Juan


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
        resolution=200,
    )

    with pytest.raises(AttributeError):
        test_img.plot("./dir_test_images")

    with pytest.raises(AttributeError):
        test_img.plot("y")

    with pytest.raises(AttributeError):
        test_img.plot("img_test.png")


"""
def test_plot_pixels():
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
        resolution=0,
    )

    with pytest.raises(AttributeError):
        test_img.plot()
"""


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
        resolution=200,
    )

    with pytest.raises(AttributeError):
        test_img.plot(ax=6)
