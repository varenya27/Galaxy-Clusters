from scipy import constants
from masses import M_total,M_gas,M_stellar
import numpy as np
from matplotlib import pyplot as plt
import csv

k = constants.Boltzmann
mu=0.61
mp = constants.proton_mass
G = constants.G
n_iter=1
val=[]
with open('tables/table-1-without_errors.csv','r') as f:
    read = csv.reader(f)
    key = next(read)
    for row in read:
        val.append(row)
plt.figure()
plt.grid()
plt.title('All Clusters')
plt.xlabel('r (kpc)')
# plt.ylabel('$M_{tot}/M_{bar} (r)$')
plt.ylabel('$M_{dark}/M_{bar} (r)$')

for cluster in val:
    name=cluster[0]
    T=float(cluster[2]) * 11.6*(10**6) #kEv to K
    beta = float(cluster[3])
    r_c = float(cluster[4])
    n_c = float(cluster[5])*(10**4)


    y=[]
    z=[]
    for i in range(0,n_iter*int(r_c)):
        r=i/n_iter
        m_tot = M_total(r,T,beta,r_c)
        m_gas = M_gas(r_c,n_c,beta,r)
        m_stellar = M_stellar(m_tot)
        m_bar = m_gas+m_stellar
        m_dark = m_tot-m_bar
        # print(m_tot/m_bar)
        y.append(m_dark/m_bar)
        z.append(m_tot/m_bar)
        
    for a,b in zip(y,z):
        print(b-a)
    
    
    x=np.arange(0,r_c)
    plt.plot(x,y)

    # plt.scatter(x, y, s=0.2,color='black')
    plt.axhline(y=1,color='teal',linestyle='-')

    # plt.savefig('figures/'+name)
plt.show()