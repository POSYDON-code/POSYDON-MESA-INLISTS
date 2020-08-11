#Y = Y_BigBang + Z/Zsolar (Ysolar - Y_BigBang)

Zsolar = 0.0142
Ysolar = 0.2703 #both from Asplund et al. 2009
#Z = 0.5 * Zsolar
Z = 1.0 * Zsolar
Y_BigBang = 0.249 #(Plank et al. 2015)

Y = Y_BigBang + Z/Zsolar * (Ysolar - Y_BigBang)

print (Z, Y)
