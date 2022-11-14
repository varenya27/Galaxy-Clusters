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
n_iter=1
val=[]
with open('tables/NCCC_all_no-errors.csv','r') as f:
    read = csv.reader(f)
    key = next(read)
    for row in read:
        val.append(row)
# plt.figure()
# plt.grid()
# # plt.title('Cluster: {}'.format(name))
# plt.title('All clusters')
# plt.xlabel('$log(a_{tot})$')
# plt.ylabel('$M_{tot}/M_{bar} (a_{tot})$')
for cluster in val:
    name=cluster[0]
    T=float(cluster[-2]) * 11.6*(10**6) #kEv to K
    beta = float(cluster[3])
    r_c = float(cluster[4])
    n_c = float(cluster[5])*(10**4)
    r_500 = float(cluster[-1])*1000

  
    y=[]
    x=[]
    y_log=[]
    x_log=[]
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
    x_axis=np.arange(start,r_500)

    for i in range(start,n_iter*int(r_500)):
        r=i/n_iter
        m_gas = M_gas(r_c,n_c,beta,r)
        m_stellar = M_stellar(M_tot,r_500,r)
        m_bar = m_gas+m_stellar
        y.append(a_tot(beta, T, r_c, r))
        y_log.append(np.log10(a_tot(beta, T, r_c, r)))
        x.append(a_bar(m_bar,r))
        x_log.append(np.log10(a_bar(m_bar,r)))
    # plt.figure()
    # plt.scatter(x,y,s=0.0001)
    plt.figure()
    plt.grid()
    plt.title('Cluster: {}'.format(name))
    # plt.title('All clusters')
    plt.xlabel('$log(a_{bar})$')
    plt.ylabel('$log(a_{tot})$')
    plt.plot(x_log,y_log)
    # plt.axhline(y = 1,color='black',linestyle = '-')
    plt.show()
    # plt.savefig('figures2/a_tot/'+name)
# plt.show()