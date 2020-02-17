import matplotlib.pyplot as plt
import numpy as np

#graph of the geodesic depending on R_s, the Schwarzschild rayon

R_s = 1


def Schwarzschild_1(r, r_0):
    return r_0/R_s - np.log( (r-R_s)/(r_0 - R_s))
   

geodesic_schwarzschild = np.vectorize(Schwarzschild_1)

R = np.linspace(0, 5, 100)
T = geodesic_schwarzschild(R, 3)


plt.figure()
plt.plot(R, T)
plt.show()