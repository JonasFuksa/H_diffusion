# Program to test the rates program

import rates
import matplotlib.pyplot as plt
import numpy as np

a = 10E-10 # lattice constant
ef = 7000 # Fermi energy in meV
dk = np.array([0, 0]) # set momentum exchange
lmbd = 5.5E-11 # Screening length
a0 = 5E-11 # Bohr radius in m
T = 300 # Temperature in K
n = 10 # number of points

e = -np.linspace(2.5, 5, n)
gamma = np.zeros(n)
f = open('transition_rates_data.txt', '+w')
for i in range(n):
    gamma[i] = rates.rate(e[i], dk, lmbd, a0, 30, ef, T)
    f.write(str(gamma[i]) + ' ')
f.close()
plt.plot(e, gamma)
plt.title('Transition rates for dk = [0,0]')
plt.xlabel('de/meV')
plt.ylabel('relative rate')
plt.savefig('transition_rates.png')
plt.show()
