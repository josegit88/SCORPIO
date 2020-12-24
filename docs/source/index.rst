
..
   Automatic created file. Don't edit



.. image:: https://travis-ci.com/josegit88/SCORPIO.svg?branch=master
   :target: https://travis-ci.com/josegit88/SCORPIO
   :alt: Build Status


.. image:: https://readthedocs.org/projects/scorpio-rdd/badge/?version=latest
   :target: https://scorpio-rdd.readthedocs.io/en/latest/?badge=latest
   :alt: Build Status


.. image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License: MIT


.. image:: https://img.shields.io/badge/python-3.7+-blue.svg
   :target: https://www.python.org/downloads/release/python-370/
   :alt: Python 3.7+


.. image:: https://img.shields.io/badge/DiSoftCompCi-FAMAF-ffda00
   :target: https://github.com/leliel12/diseno_sci_sfw
   :alt: https://github.com/leliel12/diseno_sci_sfw



.. image:: https://github.com/josegit88/SCORPIO/raw/master/res/scorpio.png
   :target: https://github.com/josegit88/SCORPIO/raw/master/res/scorpio.png
   :alt: quick tool to generate images of astrophysical objects


SCORPIO
=======

..

   Sky COllector of galaxy Pairs and Image Output


Is a tool to quick generate images of galaxy pairs, using data from different surveys.

Motivation
----------

The interacting galaxies are one of the most interesting events in the area of extragalactic astronomy. The degree of interaction can sometimes be biased by projection effects. From the osbervation point of view, a visual image obtained from the photometry of these objects makes it possible to establish and quantify some features of this interaction, such as tidal tails and other deformations caused by gravitational force.

Features
--------

The essential input parameters are AR, DEC, redshift (and other optional data such as survey and filters) for each of the galaxies of interest.

Based on these, the application has several functions that:


* Check and download the .fits files in each survey and filter, using astroquery tools.
* stack and generate a single array with the information
* Calculate the distance in physical units and in pixels.
* Export the image in the desired format and directory.

Requirements
------------

To use SCORPIO you need python 3.7 or higher, and have the following libraries:


* numpy
* astropy
* matplotlib
* astroquery:
  The latest version of astroquery can be conda or pip installed:
  conda install -c astropy astroquery
  or pip install astroquery

(It is suggested to create a virtual environment if you do not want some of these packages to be permanently hosted in your most used libraries)

Installation
------------

.. code-block::

   $ pip install scorpio-gp


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api/modules
   installation.rst
   license.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

