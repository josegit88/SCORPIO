import scorpio

#datos de prueba con dos galaxias:
[RA1, DEC1, Z1, RA2, DEC2, Z2] = [126.39162693999999, 47.296980665521900, 0.12573827000000001, 126.38991429000001, 47.305200665521902, 0.12554201000000001]

prueba_img = scorpio.gpair(ra1=RA1, dec1=DEC1, ra2=RA2, dec2=DEC2, z1=Z1, z2=Z2)

prueba_img.plot()

