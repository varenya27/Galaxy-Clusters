
import csv
from scipy import constants
import numpy as np
from scipy.integrate import quad
from astropy.constants import M_sun
from astropy.units import kg

def a_tot(beta, T,r_c,r):
    k = constants.Boltzmann
    m = constants.proton_mass
    a = (3*beta*k*T*r)/(m*(r**2+r_c**2)* (3.086*(10**19)))
    return a

def a_bar(M_bar,r):
    G = constants.G
    a = G*M_bar*(M_sun/kg)/((r* (3.086*(10**19)))**2)
    return a