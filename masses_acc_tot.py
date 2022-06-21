import csv
from scipy import constants
import numpy as np
from astropy.constants import M_sun
from astropy.units import kg
from scipy.integrate import quad
k = constants.Boltzmann
mu=0.61
mp = constants.proton_mass
G = constants.G
def M_total_acc_tot(a, Th,beta,r_c):

    const = 3*beta*k*Th/mp
    r = ((const/a)+np.sqrt((const/a)**2-4*((r_c* (3.086*(10**19)))**2)))/2
    M_r=((3*beta*k*Th*r* (3.086*(10**19)))/(G*mu*mp))*( (r/r_c)**2 / (1+(r/r_c)**2))
    return M_r*kg/M_sun

def M_gas_acc_tot(Th,r_c,n_c,beta,acc):
    def M(a):
        const = 3*beta*k*Th/mp
        rx = ((const/a)+np.sqrt((const/a)**2-4*((r_c* (3.086*(10**19)))**2)))/2
        print('*****')
        print((const/a)**2-4*r_c**2)
        n_r=n_c*((1+(rx/r_c)**2))**(-3*beta/2)*constants.atomic_mass*0.6
        return (n_r*4*np.pi*(rx* (3.086*(10**19)))**2)*(   ( const/(2*(a**2)) )* ( 1+(const/a)/np.sqrt(const**2/a**2 - 4*(r_c*((r_c* (3.086*(10**19)))**2))**2) )      )
    M_cur,err = quad(M,0,acc)
    return M_cur*kg* (3.086*(10**19))/M_sun

def M_stellar_acc_tot(m_tot):
    return (4*(10**12)*(((m_tot)/(5.7*(10**13)))**0.6))
