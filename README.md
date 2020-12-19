[![Build Status](https://travis-ci.com/josegit88/SCORPIO.svg?branch=master)](https://travis-ci.com/josegit88/SCORPIO)
[![Build Status](https://readthedocs.org/projects/scorpio-rdd/badge/?version=latest)](https://scorpio-rdd.readthedocs.io/en/latest/?badge=latest)
[![https://github.com/leliel12/diseno_sci_sfw](https://img.shields.io/badge/DiSoftCompCi-FAMAF-ffda00)](https://github.com/leliel12/diseno_sci_sfw)

![quick tool to generate images of astrophysical objects](https://raw.githubusercontent.com/josegit88/SCORPIO/master/docs/source/_static/scorpio_logo.png)

# SCORPIO (Sky COllector of galaxy Pairs and Image Output)
Is a tool to quick generate images of galaxy pairs, using data from different surveys.

## Motivation
The interacting galaxies are one of the most interesting events in the area of extragalactic astronomy. The degree of interaction can sometimes be biased by projection effects. From the osbervation point of view, a visual image obtained from the photometry of these objects makes it possible to establish and quantify some features of this interaction, such as tidal tails and other deformations caused by gravitational force.

## Features
The essential input parameters are AR, DEC, redshift (and other optional data such as survey and filters) for each of the galaxies of interest.

Based on these, the application has several functions that:
- Check and download the .fits files in each survey and filter, using astroquery tools.
- stack and generate a single array with the information
- Calculate the distance in physical units and in pixels.
- Export the image in the desired format and directory.

## Requirements
To use SCORPIO you need python 3.7 or higher, and have the following libraries:
- numpy
- astropy
- matplotlib
- astroquery:
  The latest version of astroquery can be conda or pip installed:
  conda install -c astropy astroquery
  or pip install astroquery

(It is suggested to create a virtual environment if you do not want some of these packages to be permanently hosted in your most used libraries)


## Installation
Clone this repo and then inside the local directory execute

    $ pip install -i https://test.pypi.org/simple/ scorpio==0.0.1

![quick tool to generate images of astrophysical objects](https://raw.githubusercontent.com/josegit88/SCORPIO/master/docs/source/_static/tenor.gif)
