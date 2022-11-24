import csv
from scipy import constants
import numpy as np
# from astropy.constants import M_sun
# from astropy.units import kg
from scipy.integrate import quad
m=constants.atomic_mass
M_sun=1.989e30
norm = 3.086e19


def M_total(r, T_,Beta,R_c):
    k = constants.Boltzmann
    mu=0.61
    mp = constants.proton_mass
    G = constants.G
    T,err_T=T_
    beta,err_beta=Beta
    r_c,err_r_c=R_c
    M_r=((3*beta*k*T*r* (3.086e19))/(G*mu*mp))*( (r/r_c)**2 / (1+(r/r_c)**2))
    err_M_r = np.sqrt( (M_r*err_T/T)**2 + (M_r*err_beta/beta)**2 + (M_r*(-2*err_r_c/(r_c*(1+(r/r_c)**2))))**2 )
    
    M_r=M_r/M_sun
    err_M_r = err_M_r/M_sun
    return [M_r,err_M_r]


def M_gas(R_c,N_c,Beta,lim):
    r_c,err_r_c=R_c
    n_c,err_n_c=N_c
    beta,err_beta=Beta
    def M(rx):
        n_r=n_c*((1+(rx/(r_c*norm))**2))**(-3*beta/2)*m*0.6
        return n_r*4*np.pi*(rx)**2
    M_cur,err = quad(M,0,lim*norm)
    
    def M_n_(rx):
        n_r=err_n_c*((1+(rx/(r_c*norm))**2))**(-3*beta/2)*m*0.6
        return n_r*4*np.pi*(rx)**2
    M_n,err1 = quad(M_n_,0,lim*norm)

    def M_beta_(rx):
        n_r=n_c*((1+(rx/(r_c*norm))**2))**(-3*beta/2)*m*0.6
        return n_r*4*np.pi*(rx)**2*(-1.5*np.log(1+(rx/(r_c*norm))**2))*err_beta
    M_beta,err2 = quad(M_beta_,0,lim*norm)

    def M_rc_(rx):
        n_r=n_c*((1+(rx/(r_c*norm))**2))**(-3*beta/2-1)*m*0.6*3*beta*rx**2/((r_c*norm)**3)*(err_r_c*norm)
        return n_r*4*np.pi*(rx)**2
    M_rc,err3 = quad(M_rc_,0,lim*norm)
    # print("{:e} {:e} {:e} {:e} ".format(M_cur, M_beta,M_n,M_rc))
    M_cur/=M_sun
    err_M_cur = np.sqrt((M_n)**2+(M_beta)**2+(M_rc)**2)/M_sun

    return [M_cur,err_M_cur]


def M_stellar(M_gas):
    m_gas,err_m_gas = M_gas
    if(m_gas==0): return 0,0
    # print(M_gas)
    a,err_a = 0.6,0.1
    m_star=  (4e12*(((m_gas)/(5.7e13))**a))
    dMsdMg = (4e12*(((m_gas)/(5.7e13))**a))*a/m_gas
    dMsda =  (4e12*(((m_gas)/(5.7e13))**a))*np.log(((m_gas)/(5.7e13)))
    m_star_err = np.sqrt( (dMsdMg*err_m_gas)**2 + (dMsda*err_a)**2 )
    # print("{:e} {:e} {:e}".format(m_star,dMsda*err_a,m_gas),dMsda*err_a/m_star)

    # print(dMsdMg*err_m_gas,dMsda*err_a,m_star,m_gas)
    return [m_star, m_star_err]

def M_stellar2(M_gas,R500,r):
    m_gas,err_m_gas = M_gas
    r500,err_r500 = R500
    a,err_a = 0.6,0.1
    r*=norm
    m_star = (4e12*(((m_gas)/(5.7e13))**a))*(r/r500)
    dMsdMg = (4e12*(((m_gas)/(5.7e13))**a))*(r/r500) *a/m_gas
    dMsda  = (4e12*(((m_gas)/(5.7e13))**a))*(r/r500) *np.log(((m_gas)/(5.7e13)))
    dMsdr500 = (4e12*(((m_gas)/(5.7e13))**a))*(r/(r500**2))
    # print(dMsdr500)
    
    m_star_err = np.sqrt( (dMsdMg*err_m_gas)**2 + (dMsda*err_a)**2 + (dMsdr500*err_r500)**2 )
    return [m_star, m_star_err]
