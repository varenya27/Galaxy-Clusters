from scipy import constants
from masses_acc_tot import M_total_acc_tot,M_gas_acc_tot,M_stellar_acc_tot

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

    # print(3*beta*k*T/(2*mp*r_c* (3.086*(10**19))))
    y=[]
    for i in range(1,100):

        a=i/100*(3*beta*k*T/(2*mp*r_c))
        m_tot = M_total_acc_tot(a,T,beta,r_c)
        m_gas = M_gas_acc_tot(T,r_c,n_c,beta,a)
        m_stellar = M_stellar_acc_tot(m_tot)
        m_bar = m_gas+m_stellar
        # print(m_tot/m_bar)
        y.append(m_tot/m_bar)

    print(y)
    break
    # plt.figure()
    # x=np.arange(1,100)
    # plt.plot(x,y)
    # plt.plot(x,x*0+1)
    # plt.grid()
    # plt.title('Cluster: {}'.format(name))
    # plt.xlabel('a_tot ()')
    # plt.ylabel('$M_{tot}/M_{bar} (a_tot)$')
    # plt.show()
    # break
    # plt.savefig('figures/total_acceleration'+name)