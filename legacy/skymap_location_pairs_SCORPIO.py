# -*- coding: utf-8 -*-
# Jose Benavides (jose.astroph@gmail.com)
"""
Documentation:
This program takes RA, Dec and redshift data from a sample of galaxies
or other objects in the sky. The plotSkyLoc function receives these arguments
and export a figure with the location of the data depending on the type of
selected projection and assigning a color by redshift.
"""
import numpy as np
import matplotlib.pyplot as plt
plt.close()

data_pares = np.genfromtxt("data_pares_galaxias.dat")


def plotSkyLoc(RA, Dec, redshift, org=0, title=None, projection=None,
               colormap="rainbow"):
    """RA, Dec are arrays of the same length.
    RA takes values in [0,360), Dec in [-90,90], which represent angles in
    degrees. org is the origin of the plot, 0 or a multiple of 30 degrees
    in [0,360). title is the title of the figure. projection is the kind
    of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'. colormap
    assings the color code for redshift colorbar. We can select: "rainbow",
    viridis, plasma, inferno, magma.
    """
    projection = 'mollweide' if projection is None else projection
    x = np.remainder(RA+360.-org, 360.)  # shift RA values
    ind = x > 180.
    x[ind] -= 360.  # scale conversion to [-180, 180]
    x = -x  # reverse the scale: East to the left
    tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tick_labels = np.remainder(tick_labels+360+org, 360)
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection=projection)
    color_value = redshift

    if colormap == "rainbow":
        cmap_select = plt.cm.rainbow
    if colormap == "viridis":
        cmap_select = plt.cm.viridis
    if colormap == "plasma":
        cmap_select = plt.cm.plasma
    if colormap == "inferno":
        cmap_select = plt.cm.inferno
    if colormap == "magma":
        cmap_select = plt.cm.magma

    im = ax.scatter(
        np.radians(x), np.radians(Dec), c=color_value, lw=0.4,
        marker=".", edgecolors="k", cmap=cmap_select)

    cb = fig.colorbar(
        im, ax=ax, orientation="horizontal", pad=0.08, shrink=0.75)
    cb.ax.tick_params(
        which="major", length=6, right=True, direction="in", labelsize=10)
    cb.set_label(r"redshift", fontsize=14)
    im.set_clim(0, 0.2)
    ax.set_xticklabels(tick_labels)  # we add the scale on the x axis
    if title == None:
        title = ""
    ax.set_title(title)
    ax.title.set_fontsize(11)
    ax.set_xlabel("RA")
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel("Dec")
    ax.yaxis.label.set_fontsize(12)
    ax.grid(True)

    return ax


cc = "magma"
plotSkyLoc(
    data_pares[:2001:2, 1], data_pares[:2001:2, 2], data_pares[:2001:2, 3],
    title='Sky projection of Galaxy pairs',org=0., projection='mollweide',
    colormap=cc)
plt.savefig("projection_sky_sample_y_cbar_horizontal_"+str(cc)+"_mollweide.png",
            bbox_inches='tight', dpi=200)
plt.show()
