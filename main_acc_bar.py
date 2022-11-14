from scipy import constants
from masses import M_total,M_gas,M_stellar
from accelerations import a_bar,a_tot
import numpy as np
from matplotlib import pyplot as plt
import csv

n_iter=1
k = constants.Boltzmann
mu=0.61
mp = constants.proton_mass
G = constants.G
val=[]
with open('tables/NCCC_all_no-errors.csv','r') as f:
    read = csv.reader(f)
    key = next(read)
    for row in read:
        val.append(row)
# plt.figure()
# plt.grid()
# plt.xlabel('$log(a_{bar})$')
# plt.ylabel('$M_{tot}/M_{bar} (a_{bar})$')
# plt.title('All clusters')

for cluster in val:
    name=cluster[0]
    T=float(cluster[-2]) * 11.6*(10**6) #kEv to K
    beta = float(cluster[3])
    r_c = float(cluster[4])
    n_c = float(cluster[5])*(10**4)
    r_500 = float(cluster[-1])*1000


    y=[]
    x=[]
    x_log=[]
    start=0
    M_tot = M_total(r_500,T,beta,r_c)
    for i in range(0,n_iter*int(r_500)):
        r=i/n_iter
        m_tot = M_total(r,T,beta,r_c)
        m_gas = M_gas(r_c,n_c,beta,r)
        m_stellar = M_stellar(M_tot,r_500,r)
        m_bar = m_gas+m_stellar
        if(m_tot/m_bar>=1):
            start = int(r);
            break

    for i in range(start,n_iter*int(r_500)):
        r=i/n_iter
        m_tot = M_total(r,T,beta,r_c)
        m_gas = M_gas(r_c,n_c,beta,r)
        m_stellar = M_stellar(M_tot,r_500,r)
        m_bar = m_gas+m_stellar
        m_dark = m_tot-m_bar
        y.append(m_tot/m_bar)
        # x.append(a_bar(m_bar,r))
        x_log.append(np.log10(a_bar(m_bar,r)))
        # print(x)
        # print(x,y)
    # plt.figure()
    plt.figure()
    plt.grid()
    plt.xlabel('$log(a_{bar})$')
    plt.ylabel('$M_{tot}/M_{bar} (a_{bar})$')
    plt.plot(x_log,y)
    plt.axhline(y = 1,color='black',linestyle = '-')
    plt.title('Cluster: {}'.format(name))
    plt.savefig('figures2/a_bar/'+name)
    # plt.savefig('figures/all/a_bar_logarithmic')
    plt.show()
    # break