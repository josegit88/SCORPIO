# -*- coding: utf-8 -*-
# Jose Benavides (jose.astroph@gmail.com)

# This file is part of the
#   SCORPIO Project (https://github.com/josegit88/SCORPIO).
# Copyright (c) 2020, Jose Benavides
# License: MIT
#   Full Text: https://github.com/josegit88/SCORPIO/blob/add-license-2/LICENSE

import scorpio
import matplotlib.pyplot as plt

# data of both galaxies:
[RA1, DEC1, Z1, RA2, DEC2, Z2] = [
    234.47982166000000,
    27.915027615521904,
    0.13499665999999999,
    234.48308671000001,
    27.913794015521905,
    0.13481650000000001,
]

SURV = "SDSS"
prueba_img = scorpio.gpair(
    ra1=RA1,
    dec1=DEC1,
    ra2=RA2,
    dec2=DEC2,
    z1=Z1,
    z2=Z2,
    survey=SURV,
    resolution=500,
)

prueba_img.plot(
    save_Img=True,
    imgName="pair_gp_" + str(SURV) + ".png",
    size=10,
    color_map="bone",
)
