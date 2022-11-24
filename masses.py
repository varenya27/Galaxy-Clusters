import csv
from scipy import constants
import numpy as np
from astropy.constants import M_sun
from astropy.units import kg
from scipy.integrate import quad

def M_total(r, Th,beta,r_c):
    k = constants.Boltzmann
    mu=0.61
    mp = constants.proton_mass
    G = constants.G
    M_r=((3*beta*k*Th*r* (3.086*(10**19)))/(G*mu*mp))*( (r/r_c)**2 / (1+(r/r_c)**2))
    return M_r*kg/M_sun

def M_gas(r_c,n_c,beta,lim):
    def M(rx):
        n_r=n_c*((1+(rx/r_c)**2))**(-3*beta/2)*constants.atomic_mass*0.6
        return n_r*4*np.pi*(rx* (3.086*(10**19)))**2
    M_cur,err = quad(M,0,lim)
    return M_cur*kg* (3.086*(10**19))*2/M_sun

def M_stellar(m_gas,r_500,r):
    return (r/r_500)*(4*(10**12)*(((m_gas)/(5.7*(10**13)))**0.6))

def M_gas1(r_c,n_c,beta,lim):
    norm = 3.086e19
    m=constants.atomic_mass
    # r_c,err_r_c=R_c
    # n_c,err_n_c=N_c
    # beta,err_beta=Beta
    def M(rx):
        n_r=n_c*((1+(rx/(r_c*norm))**2))**(-3*beta/2)*m*0.6
        return n_r*4*np.pi*(rx)**2
    M_cur,err = quad(M,0,lim*norm)
    return M_cur*kg/M_sun