from scipy import constants
from masses import M_total,M_gas,M_stellar
from accelerations import a_bar,a_tot
import numpy as np
from matplotlib import pyplot as plt
import csv

k = constants.Boltzmann
mu=0.61
mp = constants.proton_mass
G = constants.G
val=[]
with open('table-1-without_errors.csv','r') as f:
    read = csv.reader(f)
    key = next(read)
    for row in read:
        val.append(row)

for cluster in val:
    name=cluster[0]
    T=float(cluster[2]) * 11.6*(10**6) #kEv to K
    beta = float(cluster[3])
    r_c = float(cluster[4])
    n_c = float(cluster[5])*(10**4)


    y=[]
    x=[]
    for i in range(0,100*int(r_c)):
        r=i/100
        m_tot = M_total(r,T,beta,r_c)
        m_gas = M_gas(r_c,n_c,beta,r)
        m_stellar = M_stellar(m_tot)
        m_bar = m_gas+m_stellar
        # print(m_tot/m_bar)
        y.append(m_tot/m_bar)
        x.append(a_bar(m_bar,r))
        # print(x)
        # print(x,y)
    plt.figure()
    x1=[1]*len(x)
    plt.plot(x,y)
    plt.axhline(y = 1,color='orange',linestyle = '-')
    plt.grid()
    plt.title('Cluster: {}'.format(name))
    plt.xlabel('$a_{bar}$ ($ms^{-2}$)')
    plt.ylabel('$M_{tot}/M_{bar} (a_{bar})$')
    plt.savefig('figures/baryonic_acceleration/'+name)
    # plt.show()
    # break