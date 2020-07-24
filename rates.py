# This is a program to calculate transition rates of a hydrogen atom on the surface of a metal

# This version takes into account only the bulk electrons of the metal, treated as a Fermi gas at a fixed temperature

# To compute the scattering cross section, first Born approximation is used for an electron scattering of a screened Hydrogen atom in its ground state (right now the ground state of the screened atom is assumed to have the same functional form as in the full Coulomb case, but with variable Bohr radius)

# Energy unit is meV

# Maximum band index should be adjusted depending on the Fermi energy of the metal considered and the temperature

import numpy as np

h = 6.58E-13 # Planck const. in meV
kb = 8.62E-2 # Boltzmann constant in meV
m = 9.11E-31 # Mass of the electron
omega = 10000.0/h # Frequency of the harmonic oscillator, with which the z-potential is calculated
lz_sq = 4.0*h/(m*omega) # Length scale of the harmonic oscillator

# This is the main rates function. 
# Imputs: change of hydrogen energy de, change of hydrogen k vector dk, screening length lmbd, Bohr radius of the screened atom a0, number of states considered in one direction in the 1st Brillouin zone d, Fermi energy ef, temperature T
# Outputs: relative transition rate for this transition
def rate(de, dk, lmbd, a0, d, ef, T):
    # Set up k-space
    k = np.zeros((2*d, 2*d, 2*d, 3))
    # Limits on initial state imposed by Fermi-Dirac statistics
    e_min = ef - kb*T + de
    e_max = ef + kb*T
    cross_section = 0.0
    if e_min < e_max and e_min > 0:
        k_min = np.sqrt(2*m*e_min)/h
        k_max = np.sqrt(2*m*e_max)/h
        k1d = np.concatenate((np.linspace(k_min/np.sqrt(5), k_max, d), np.linspace (-k_max, -k_min/np.sqrt(5), d)))
        # all possible initial states
        for i in range(2*d):
            for j in range(2*d):
                for l in range(2*d):
                    k[i,j,l,0] = k1d[i]
                    k[i,j,l,1] = k1d[j]
                    k[i,j,l,2] = k1d[l]
        for i in range(2*d):
            for j in range(2*d):
                for l in range(2*d):
                    k_init = k[i,j,l]
                    e_init = e(k_init)
                    e_fin = e_init - de
                    if e_fin > 0 and e_init > e_min and e_init < e_max:
                        # Calculate the change in kz to conserve energy and x-y momentum
                        mod_k_fin = np.sqrt(2*m*e_fin)/h
                        kz_sq = mod_k_fin**2 - (k_init[0] - dk[0])**2 - (k_init[1] - dk[1])**2
                        if kz_sq > 0:
                            q1 = np.array([dk[0], dk[1], k_init[2] + np.sqrt(kz_sq)])
                            q2 = np.array([dk[0], dk[1], k_init[2] - np.sqrt(kz_sq)])
                            cross_section += n(e_init, ef, T) * (1 - n(e_fin, ef, T)) * crssct(q1, lmbd, a0)
                            cross_section += n(e_init, ef, T) * (1 - n(e_fin, ef, T)) * crssct(q2, lmbd, a0)
    print(cross_section)
    return cross_section
                    



# Function to evaluate the Fermi-dirac occupancy of a state of energy e, given Fermi energy (chamical potential) ef and temperature T
def n(e, ef, T):
    return 1/(1+np.exp((e-ef)/kb*T))

# Function to calculate the relative cross section of a transition with change of k vector q, given screening length lmbd and Bohr radius a0
def crssct(q0, lmbd, a0):
    q = np.sum(np.multiply(q0, q0))
    sct = 0.0
    if q != 0:
        sct = ((16*np.pi/(q**2*a0*((2/a0)**2 + q**2)**2) - 1/(1/lmbd**2 + q**2))*np.exp(-q0[2]**2*lz_sq/16))**2
    return sct
# Dispersion of bulk electrons; returns energy of an electron of wavevector k
# At the moment free electrons, but can be changed later for more realistic situations
def e(k):
    return h**2 * np.sum(np.multiply(k, k))/(2*m)
