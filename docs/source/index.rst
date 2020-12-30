
..
   Automatic created file. Don't edit


SCORPIO
=======

Sky COllector of galaxy Pairs and Image Output (Scorpio) Is a tool to quick analyze and generate images of galaxy pairs, using data from different surveys.


.. image:: https://github.com/josegit88/SCORPIO/raw/master/res/scorpio.png
   :target: https://github.com/josegit88/SCORPIO/raw/master/res/scorpio.png
   :alt: logo-scorpio


----


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


Motivation
----------

The interacting galaxies are one of the most interesting events in the area of extragalactic astronomy. The degree of interaction can sometimes be biased by projection effects. From the osbervation point of view, a visual image obtained from the photometry of these objects makes it possible to establish and quantify some features of this interaction, such as tidal tails and other deformations caused by gravitational force.

Features
--------

The essential input parameters are RA, DEC, redshift (and other optional data such as survey and filters) for each of the galaxies of interest.

Based on these, the application has several functions that:


* Check and download the ``.fits`` files in each survey and filter, using astroquery tools.
* stack and generate a single array with the information
* Calculate the distance in physical units and in pixels.
* Export the image in the desired format and directory.

Code and issues
---------------

The entire source code of is hosted in GitHub
`https://github.com/josegit88/SCORPIO <https://github.com/josegit88/SCORPIO>`_

License
-------

SCORPIO is under
`MIT <https://www.tldrlegal.com/l/mit>`_

A short, permissive software license. Basically, you can do whatever you want as long as you include the original copyright and license notice in any copy of the software/source.  There are many variations of this license in use.

Installation
------------

This is the recommended way to install SCORPIO.

Installing  with pip
^^^^^^^^^^^^^^^^^^^^

Make sure that the Python interpreter can load SCORPIO code.
The most convenient way to do this is to use virtualenv, virtualenvwrapper, and pip.

After setting up and activating the virtualenv, run the following command:

.. code-block:: console

   $ pip install scorpio-gp
   ...

That should be it all.

Installing the development version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you’d like to be able to update your SCORPIO code occasionally with the latest bug fixes and improvements, follow these instructions:

Make sure that you have Git installed and that you can run its commands from a shell.
(Enter *git help* at a shell prompt to test this.)

Check out SCORPIO main development branch like so:

.. code-block:: console

   $ git clone https://github.com/josegit88/SCORPIO.git
   ...

This will create a directory *SCORPIO* in your current directory.

Then you can proceed to install with the commands

.. code-block:: console

   $ cd SCORPIO
   $ pip install -e .
   ...

Documentation
-------------

The full documentation of the project are available in
`https://scorpio-rdd.readthedocs.io <https://scorpio-rdd.readthedocs.io/>`_

Contact
-------

For bugs or question please contact

..

   **José Benavides** `jose.astroph@gmail.com <jose.astroph@gmail.com>`_



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   tutorial.ipynb
   api.rst



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

