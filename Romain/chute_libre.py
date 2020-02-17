import sympy as s
import math as m
import numpy as np
import matplotlib.pyplot as plt
from sympy import init_printing
import scipy.integrate as sc


c = 1
R_s = 1

def F_diff(Y, t, r_0 = 0.5):
    return np.array([Y[1], -c**2 /2 /(Y[0]**2), R_s/c * np.sqrt(1-1/r_0) / (1-1/Y[0])])

time = np.linspace(0, 25, 30000)

Y = sc.odeint(F_diff, np.array([5, 0, 0]), time)

plt.figure()
plt.subplot(1,2,1)
plt.grid()
plt.title("Chute libre d'une particule vers un trou noir de Schwarschild, R_s = 5")
plt.ylabel("Rayon R/R_s")
plt.xlabel("Temps propre tau/(R_S /c)")
plt.plot(time, Y[:,0], 'r')
plt.subplot(1,2,2)
plt.grid()
plt.plot(Y[:,2], Y[:,0], 'b', label = "integration du temps propre")
plt.legend()
plt.show()