import csv
from scipy import constants
import numpy as np
# from astropy.constants import M_sun
# from astropy.units import kg
from scipy.integrate import quad

M_sun=1.989e30
def M_total(r, Th,beta,r_c):
    k = constants.Boltzmann
    mu=0.61
    mp = constants.proton_mass
    G = constants.G
    Th,err_Th=Th
    beta,err_beta=beta
    r_c,err_r_c=r_c
    M_r=((3*beta*k*Th*r* (3.086e19))/(G*mu*mp))*( (r/r_c)**2 / (1+(r/r_c)**2))
    err_M_r = M_r*err_Th/Th + M_r*err_beta/beta + M_r*(-2*err_r_c/(r_c*(1+(r/r_c)**2)))
    
    M_r=M_r/M_sun
    err_M_r = err_M_r/M_sun
    return [M_r,err_M_r]

def M_gas(R_c,N_c,Beta,lim):
    r_c,err_r_c=R_c
    r_c,err_r_c,lim=r_c*3.086e19,err_r_c*3.086e19,lim*3.086e19
    n_c,err_n_c=N_c
    beta,err_beta=Beta
    def M(rx):
        n_r=n_c*((1+(rx/r_c)**2))**(-3*beta/2)*constants.atomic_mass*0.6
        return n_r*4*np.pi*(rx)**2
    M_cur,err = quad(M,0,lim)

    def M_err1(rx):
        n_r=err_n_c*((1+(rx/r_c)**2))**(-3*beta/2)*constants.atomic_mass*0.6
        return n_r*4*np.pi*(rx)**2
    M_n1,err1 = quad(M,0,lim)


    def M_err2(rx):
        n_r=n_c*((1+(rx/r_c)**2))**(-3*beta/2)*constants.atomic_mass*0.6
        return n_r*4*np.pi*(rx)**2*(-1.5*np.log(1+(r/r_c)**2))
    M_beta,err2 = quad(M,0,lim)
    M_beta*= err_beta
    err2*= err_beta

    def M_err3(rx):
        n_r=n_c*((1+(rx/r_c)**2))**(-3*beta/2)*constants.atomic_mass*0.6
        return n_r*4*np.pi*(rx)**2*(3*beta)*((rx/r_c)**2)/((rx/r_c)+1)**(3*beta/2+1)
    M_r_c,err3 = quad(M,0,lim)
    M_r_c*= err_r_c/r_c
    err2*= err_r_c/r_c

    # M_rc = 4*np.pi*0.6*constants.atomic_mass*n_c*(r_c**2)*(2**(-1.5*beta))*err_r_c*(3.086e19)**3
    M_cur*=1/M_sun
    err_M_cur = (M_n1+M_beta+M_r_c+err1+err2+err3)/M_sun
    print("{:e} {:e}".format(M_cur,err_M_cur))
    return [M_cur,err_M_cur]

def M_stellar(M_gas):
    m_gas,err_m_gas = M_gas
    if(m_gas==0): return 0,0
    # print(M_gas)
    a,err_a = 0.6,0.1
    m_star=  (4*(10**12)*(((m_gas)/(5.7*(10**13)))**a))
    dMsdMg = (4*(10**12)*(((m_gas)/(5.7*(10**13)))**a))*a/m_gas
    dMsda =  (4*(10**12)*(((m_gas)/(5.7*(10**13)))**a))*np.log(((m_gas)/(5.7*(10**13))))
    m_star_err = dMsdMg*err_m_gas+dMsda*err_a
    # print(dMsdMg*err_m_gas,dMsda*err_a,m_star,m_gas/1e13)
    return [m_star, m_star_err]