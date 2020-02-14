# -*- coding: utf-8 -*-
# Jose Benavides (jose.astroph@gmail.com)
"""
Documentation
Este programa toma los datos de RA, Dec y redshift de un conjunto de pares de
galaxias u otros objetos y calcula las distancias entre ellas y al observador.
Exporta una figura con el histograma de las distancias entre las galaxias,
otro histograma de fracción de galaxias a una deteminada distancia en escala
logarítmica y delimitando por regiones según el tipo de interacción (fuerte,
media o débil) considerada a partir de su distancia relativa. También exporta
un gráfico de la distancia entre las galaxias y la distancia al observador.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.close()

data_pares = np.genfromtxt("small_sample_data_galaxy_pairs.dat")
list_galA = range(1, len(data_pares)+1, 2)
list_galB = range(0, len(data_pares), 2)

H0 = 73.52  # Hubble constant
c_luz = 3.e5  # speed of light in km/s
D = (c_luz/H0)*data_pares[:, 3]  # estimated distance to the galaxies

max(data_pares[:, 3])
min(data_pares[:, 3])

dist_pairs = []

for ii in range(len(list_galA)):
    galA = list_galA[ii]
    galB = list_galB[ii]

    # coords of DEC y RA for each pair:
    data_pares[galA, 2], data_pares[galA, 1]
    data_pares[galB, 2], data_pares[galB, 1]

    # projected angular distance determination:
    theta_AB = np.sqrt(
        (data_pares[galA, 2] - data_pares[galB, 2])**2 +
        (data_pares[galA, 1] - data_pares[galB, 1])**2)

    # conversion to radians:
    theta_rad = theta_AB*np.pi/180.

    # distance between galaxies in kpc, taken as the arc length:
    S_AB = (np.mean((D[galA], D[galB]))*theta_rad)*1000
    dist_pairs.append((ii, theta_AB, S_AB, np.mean((D[galA], D[galB]))))

dist_pairs = np.array(dist_pairs)
np.savetxt("distancia_entre_pares.txt", dist_pairs,
           fmt=["%4.f", "%8.4f", "%8.2f", "%8.2f"])
# 0:pair number
# 1:angle between pair in degrees
# 2:distance between galaxies in kpc
# 3:distance to observer in Mpc


# =============== plots: ================
fig = plt.figure(figsize=(6, 6))

llss = 10
tMl = 5
tml = 3
plt.locator_params(axis="y", tight=True, nbins=6)
plt.tick_params(which="major", length=tMl, top=True, left=True, right=True,
                direction="in", labelsize=llss)

plt.minorticks_on()

plt.tick_params(which="minor", length=tml, top=True, left=True, right=True,
                direction="in", labelsize=llss)

plt.locator_params(axis="x", tight=True, nbins=8)
plt.locator_params(axis="y", tight=True, nbins=5)

plt.hist(dist_pairs[:, 2], histtype="step", log=False, bins=14,
         range=[0, 700], color="red", linestyle="-", hatch="/")
plt.axvline(x=np.mean(dist_pairs[:, 2]), color="b", ls="--", lw=2)
plt.xlabel(r"$R_{AB}$ [kpc]", fontsize=14)
plt.ylabel(r"N", fontsize=14)
plt.xlim(0, 700)

plt.savefig("hist_distances_between_pairs.png", bbox_inches='tight', dpi=200)
plt.show()
# plt.close()

# -------------------------------------------------
fig = plt.figure(figsize=(6, 6))

llss = 10
tMl = 5
tml = 3
plt.locator_params(axis="y", tight=True, nbins=6)
plt.tick_params(which="major", length=tMl, top=True, left=True, right=True,
                direction="in", labelsize=llss)

plt.minorticks_on()
plt.tick_params(which="minor", length=tml, top=True, left=True, right=True,
                direction="in", labelsize=llss)

plt.locator_params(axis="x", tight=True, nbins=8)
plt.locator_params(axis="y", tight=True, nbins=8)

# vertical regions for type of interaction:
media_R_AB = np.mean(np.log10(dist_pairs[:, 2]))
std_R_AB = np.std(np.log10(dist_pairs[:, 2]))
plt.axvspan(0, media_R_AB - std_R_AB, facecolor='g', alpha=0.7)  # Strong
plt.axvspan(media_R_AB - std_R_AB, media_R_AB + std_R_AB,
            facecolor='g', alpha=0.4)  # Medium
plt.axvspan(media_R_AB + std_R_AB, 3, facecolor='g', alpha=0.2)  # weak

pesos = np.ones_like(
        np.log10(dist_pairs[:, 2]))/float(len(np.log10(dist_pairs[:, 2])))
plt.hist(np.log10(dist_pairs[:, 2]), histtype="step", weights=pesos, log=False,
         bins=12, range=[0, 3], color="red", linestyle="-", hatch="/")

plt.axvline(x=np.mean(np.log10(dist_pairs[:, 2])), color="b", ls="--", lw=2)
plt.xlabel(r"$log(R_{AB})$ [kpc]", fontsize=14)
plt.ylabel(r"$f = \frac{N}{N_{tot}}$", fontsize=14)
plt.xlim(0.1, 3)

prom_SAB = float("{0:.2f}".format(10**(np.mean(np.log10(dist_pairs[:, 2])))))

plt.text(0.5, 0.32, r"$\bar{R}_{AB}$ = "+str(prom_SAB)+" kpc",
         color="b", fontsize=13)

name_hist_norm = "HistDistances_log_norm_and_interaction_Type.png"
plt.savefig(name_hist_norm, bbox_inches='tight', dpi=200)
plt.show()
# plt.close()

# =========== Distance between pair vs distance at observer: =================
fig = plt.figure(figsize=(6, 6))
gs = gridspec.GridSpec(2, 2, hspace=0.0, wspace=0.0,
                       width_ratios=[4, 1], height_ratios=[1, 4])

ax1 = plt.subplot(gs[0])  # Horizontal histogram
ax2 = plt.subplot(gs[1])  # empty
ax3 = plt.subplot(gs[2])  # main plot
ax4 = plt.subplot(gs[3])  # vertical histogram

llss = 10
tMl = 5
tml = 3
ax3.locator_params(axis="y", tight=True, nbins=6)
ax1.tick_params(which="major", length=tMl, labelbottom=False, top=True,
                left=True, right=True, direction="in", labelsize=llss)
ax3.tick_params(which="major", length=tMl, top=True, left=True, right=True,
                direction="in", labelsize=llss)
ax4.tick_params(which="major", labelleft=False, length=tMl, top=True,
                left=True, right=True, direction="in", labelsize=llss)

ax1.minorticks_on()
ax3.minorticks_on()
ax4.minorticks_on()
ax1.tick_params(which="minor", length=tml, top=True, left=True, right=True,
                direction="in", labelsize=llss)
ax3.tick_params(which="minor", length=tml, top=True, left=True, right=True,
                direction="in", labelsize=llss)
ax4.tick_params(which="minor", length=tml, top=True, left=True, right=True,
                direction="in", labelsize=llss)

ax1.locator_params(axis="y", tight=True, nbins=4)
ax4.locator_params(axis="x", tight=True, nbins=4)
ax3.locator_params(axis="x", tight=True, nbins=8)  # central plot
ax3.locator_params(axis="y", tight=True, nbins=8)  # central plot

ax3.plot(np.log10(dist_pairs[:, 2]), np.log10(dist_pairs[:, 3]), color="k",
         marker="o", fillstyle="none", ms=4, ls="none")

# upper histograms:
ax1.hist(np.log10(dist_pairs[:, 2]), histtype="step", log=False, bins=12,
         range=[0, 3], color="red", linestyle="-", hatch="/")

# lateral histograms:
ax4.hist(np.log10(dist_pairs[:, 3]), histtype="step", log=False, bins=12,
         range=[1.5, 3], color="red", linestyle="-",
         hatch="/", orientation="horizontal")

ax3.set_xlabel(r"$log(R_{AB})$ [kpc]", fontsize=14)
ax3.set_ylabel(r"$log(D_{AB})$ [Mpc]", fontsize=14)
ax1.set_ylabel("N", fontsize=14)
ax4.set_xlabel("N", fontsize=14)
ax1.axis([0, 2.7, 0, 1600])
ax3.axis([0, 2.7, 1.5, 3.05])
ax4.axis([0, 1200, 1.5, 3.05])
ax2.axis('off')
name_plot3 = "distance_between_pairs_and_observers_with_histograms.png"
plt.savefig(name_plot3, bbox_inches='tight', dpi=200)
plt.show()
# plt.close()

# ---------------------------------
print("-------------------\nThe program has ended successfully")
