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
with open('tables/NCCC_all_no-errors.csv','r') as f:
    read = csv.reader(f)
    key = next(read)
    for row in read:
        val.append(row)
# plt.figure()
# plt.grid()
# plt.title('All Clusters')
# plt.xlabel('r (kpc)')
# plt.ylabel('$M_{tot}/M_{bar} (r)$')
plt.ylabel('$M_{dark}/M_{bar} (r)$')

for cluster in val:
    name=cluster[0]
    T=float(cluster[-2]) * 11.6*(10**6) #kEv to K
    beta = float(cluster[3])
    r_c = float(cluster[4])
    n_c = float(cluster[5])*(10**4)
    r_500 = float(cluster[-1])*1000
    y,z,m1,m2,m3=[],[],[],[],[]
    
    M_tot = M_total(r_500,T,beta,r_c)
    start=0
    for i in range(0,n_iter*int(r_500)):
        r=i/n_iter
        m_tot = M_total(r,T,beta,r_c)
        m_gas = M_gas(r_c,n_c,beta,r)
        m_stellar = M_stellar(M_tot,r_500,r)
        m_bar = m_gas+m_stellar
        if(m_tot/m_bar>1):
            start = int(r);
            print(1)
            break

    for i in range(start,n_iter*int(r_500)):
        r=i/n_iter
        m_tot = M_total(r,T,beta,r_c)
        m_gas = M_gas(r_c,n_c,beta,r)
        m_stellar = M_stellar(M_tot,r_500,r)
        m_bar = m_gas+m_stellar
        m_dark = m_tot-m_bar

        # m1.append(m_tot)
        # m2.append(m_gas)
        # m3.append(m_stellar)
        y.append(m_dark/m_bar)
        z.append(m_tot/m_bar)
    # writer.writerow[name, m1[-1]/(10**14),m2[-1]/(10**13),m3[-1]/(10**13)]
    # for a,b in zip(y,z):
    #     print(b-a)
    plt.figure()
    plt.grid()
    plt.title('Cluster '+name)
    plt.xlabel('r (kpc)')
    x=np.arange(start,r_500)
    plt.plot(x,y)
    plt.plot(x,z)
    # plt.plot(x,m1)
    # plt.plot(x,m2)
    # plt.plot(x,m3)

    plt.legend(['$M_{dark}/M_{bar}$','$M_{tot}/M_{bar}$'])
    # plt.legend(['$M_{tot}$','$M_{gas}$','$M_{stellar}$'])

    # # plt.scatter(x, y, s=0.2,color='black')
    # plt.axhline(y=1,color='teal',linestyle='-')

    plt.savefig('figures2/radius/'+name)
    # plt.show()
