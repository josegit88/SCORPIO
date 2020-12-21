Installation
============


This is the recommended way to install scorpio.

Installing  with pip
^^^^^^^^^^^^^^^^^^^^

To use SCORPIO you need python 3.7 or higher, and have the following libraries:

numpy

astropy

matplotlib

astroquery: The latest version of astroquery can be conda or pip installed: conda install -c astropy astroquery or pip install astroquery

(It is suggested to create a virtual environment if you do not want some of these packages to be permanently hosted in your most used libraries)

After setting up and activating the virtualenv, run the following command:

.. code-block:: console

   $ pip install -i https://test.pypi.org/simple/ scorpio==0.0.1
   ...

That should be it all.



Installing the development version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you’d like to be able to update your SCORPIO code occasionally with the
latest bug fixes and improvements, follow these instructions:

Make sure that you have Git installed and that you can run its commands from a shell.
(Enter *git help* at a shell prompt to test this.)

Check out SCORPIO main development branch like so:

.. code-block:: console

   $ git clone https://github.com/josegit88/SCORPIO.git
   ...

This will create a directory *SCORPIO* in your current directory.


