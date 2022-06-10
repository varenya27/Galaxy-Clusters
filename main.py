from scipy import constants
from masses import M_total,M_gas,M_stellar
import numpy as np
from matplotlib import pyplot as plt

k = constants.Boltzmann
mu=0.61
mp = constants.proton_mass
G = constants.G

T=3.01* 11.6*(10**6)
beta=0.575
r_c=33
n_c=5.47*(10**4)

# r = np.linspace(0, 100000)
if(1):
    r=r_c
    m_tot = M_total(r,T,beta,r_c)
    m_gas = M_gas(r_c,n_c,beta,r)
    m_stellar = M_stellar(m_tot)
    # print(r,'{:e}'.format(m_tot),'{:e}'.format(m_gas),'{:e}'.format(m_stellar))
    # print(m_tot/(m_gas+m_stellar))
y=[]
for i in range(0,10*int(r_c)):
    r=i/10
    m_tot = M_total(r,T,beta,r_c)
    m_gas = M_gas(r_c,n_c,beta,r)
    m_stellar = M_stellar(m_tot)
    m_bar = m_gas+m_stellar
    print(m_tot/m_bar)
    y.append(m_tot/m_bar)
plt.figure()
x=np.linspace(0,1,10*int(r_c))
plt.plot(x,y)
plt.plot(x,x*0+1)
plt.grid()
plt.show()