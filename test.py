# -*- coding: utf-8 -*-
# Jose Benavides (jose.astroph@gmail.com)
"""
Test about the scorpio functions
"""

import urllib

import numpy as np

import pytest

from astroquery.skyview import SkyView
from astropy.io import fits

import scorpio



def test_download_and_stack(monkeypatch):
    aa = [126.39162693999999, 47.296980665521900, 0.12573827000000001]
    bb = [126.38991429000001, 47.305200665521902, 0.12554201000000001]

    def mock_func_g1g2(position, survey, **kwargs):
        pos = [position.ra.value, position.dec.value]

        if pos == aa[:-1] and survey == "SDSSg":
            path = [fits.open("test_data/SDSS_image_0_filter_g.fits")]
        elif pos == aa[:-1] and survey == "SDSSi":
            path = [fits.open("test_data/SDSS_image_0_filter_i.fits")]
        elif pos == bb[:-1] and survey == "SDSSg":
            path = [fits.open("test_data/SDSS_image_1_filter_g.fits")]
        elif pos == bb[:-1] and survey == "SDSSi":
            path = [fits.open("test_data/SDSS_image_1_filter_i.fits")]

        return path

    monkeypatch.setattr(SkyView, "get_images", mock_func_g1g2)

    expected_g1 = np.array([
        [-2.45361328e-02, -6.26220703e-02, -1.12152100e-02],
        [ 1.62239075e-02,  2.85937500e+00,  1.81388855e-02],
        [-3.85131836e-02,  8.28247070e-02,  8.54492188e-04]])

    expected_g2 = np.array([
        [ 1.13616943e-01, -1.66931152e-02,  4.88281250e-04],
        [ 1.61819458e-02,  2.73291016e+00,  9.80529785e-02],
        [ 2.47344971e-02, -2.17952728e-02,  4.67758179e-02]])

    stack_g1, stack_g2 = scorpio.stack_pair(aa, bb, plx=3)

    np.testing.assert_allclose(stack_g1, expected_g1, rtol = 1e-5)
    np.testing.assert_allclose(stack_g2, expected_g2, rtol = 1e-5)


def test_download_invalid_filter():
    aa = [126.39162693999999, 47.296980665521900, 0.12573827000000001]
    bb = [126.38991429000001, 47.305200665521902, 0.12554201000000001]
    with pytest.raises(ValueError):
        scorpio.stack_pair(aa, bb, plx=200, filters=["u","t"])


def test_no_filter_to_stack(monkeypatch):
    def mock_func(*args, **kwargs):
        raise urllib.error.HTTPError(
            url="http://from.mock", code=404, msg="from mock",
            hdrs=None, fp=None)

    monkeypatch.setattr(SkyView, "get_images", mock_func)

    cc = [229.38791793999997, -15.1525941155219059, 0.12589231000000001]
    dd = [229.38705890000003, -15.1513408255219060, 0.12634666999999999]
    with pytest.raises(scorpio.NoFilterToStackError):
        scorpio.stack_pair(cc, dd, plx=2450, filters=["g"])


def test_stack_code_error(monkeypatch):
    def mock_func(*args, **kwargs):
        raise urllib.error.HTTPError(
            url="http://from.mock", code=500, msg="from mock",
            hdrs=None, fp=None)

    monkeypatch.setattr(SkyView, "get_images", mock_func)

    cc = [229.38791793999997, -15.1525941155219059, 0.12589231000000001]
    dd = [229.38705890000003, -15.1513408255219060, 0.12634666999999999]
    with pytest.raises(urllib.error.HTTPError):
        scorpio.stack_pair(cc, dd, plx=2450, filters=["g"])
