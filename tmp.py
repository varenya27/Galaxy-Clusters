from scipy import constants
from astropy.units import kg
from astropy.constants import M_sun

def M_total(r, Th,beta,r_c):
    k = constants.Boltzmann
    mu=0.61
    mp = constants.proton_mass
    G = constants.G
    M_r=((3*beta*k*Th*r)/(G*mu*mp))*( (r/r_c)**2 /( 1+(r/r_c)**2))
    return M_r*kg/M_sun

T=3.01* 11.6*(10**6)
beta=0.575
r_c=33* (3.086*(10**19))
n_c=5.47*(10**4)
r=r_c
print(( (r/r_c)**2 / (1+(r/r_c)**2)))
print(M_total(r_c,T,beta,r_c))