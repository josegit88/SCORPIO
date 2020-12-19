.. SCORPIO documentation master file, created by
   sphinx-quickstart on Tue Dec  8 12:43:51 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SCORPIO's documentation!
===================================

.. Escribir una descripcion corta del proyecto aca:

**SCORPIO** (Sky COllector of galaxy Pairs and Image Output) It is an application that generates images of pairs of galaxies, in different surveys and filters.

Is a tool to quick generate images of galaxy pairs, using data from different surveys.
To use SCORPIO you need python 3.7 or higher, and have the following libraries:

numpy

astropy

matplotlib

astroquery: The latest version of astroquery can be conda or pip installed: conda install -c astropy astroquery or pip install astroquery

(It is suggested to create a virtual environment if you do not want some of these packages to be permanently hosted in your most used libraries)

| **Authors**
| Jose Benavides (E-mail: jose.benavides@gmail.com)


This program receives RA and Dec information from a pair of galaxies.
interacting (or nearby) or other individual objects and download the data
Corresponding .fits of survey data and list of filters:

******************************************************

        Survey              >>>      Filters

******************************************************

Sloan Digital sky Survey (SDSS) >>> [u, g, r, i, z]

Two Micron All Sky Survey (2MASS) >>> [J, H, K]

Wide Field Infrared Survey Explorer (WISE) >>> [3.4, 4.6, 12, 22]

Later it processes and exports an image.

For this it has a series of steps that use some specific purpose functions:

******************************************************

        Process              >>>      Function

******************************************************

Download .fits data          >>>    download_data

Stack the information        >>>      stack_pair

Calculate distance to the
observer Mpc and the
separation between the pair  >>>      distances

Generate the image           >>>      Image.plot


--------------------------------------------------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api/modules


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
