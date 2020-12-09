SCORPIO
=======

.. Escribir alguna documentacion general y corta aca.

The essential input parameters are AR, DEC, redshift (and other optional ones such as survey and filters) for each of the galaxies of interest.

Based on these, the application has several functions that:

-Check and download the .fits files in each survey and filter, using astroquery tools.

-stack and generate a single array with the information

-Calculate the distance in physical units and in pixels.

-Export the image in the desired format and directory.

.. toctree::
   :maxdepth: 4
   
   scorpio
