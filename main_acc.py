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
n_iter=10
val=[]
with open('tables/table-1-without_errors.csv','r') as f:
    read = csv.reader(f)
    key = next(read)
    for row in read:
        val.append(row)
plt.figure()
plt.grid()
# plt.title('Cluster: {}'.format(name))
plt.title('All clusters')
plt.xlabel('$log(a_{bar})$')
plt.ylabel('$log(a_{tot})$')
for cluster in val:
    name=cluster[0]
    T=float(cluster[2]) * 11.6*(10**6) #kEv to K
    beta = float(cluster[3])
    r_c = float(cluster[4])
    n_c = float(cluster[5])*(10**4)


    y=[]
    x=[]
    x_log=[]
    y_log=[]
    for i in range(1,n_iter*int(r_c)):
        r=i/n_iter
        y.append(a_tot(beta, T,r_c,r))
        y_log.append(np.log10(a_tot(beta,T,r_c,r)))
    # plt.figure()
        m_tot = M_total(r,T,beta,r_c)
        m_gas = M_gas(r_c,n_c,beta,r)
        m_stellar = M_stellar(m_tot)
        m_bar = m_gas+m_stellar
        x.append(a_bar(m_bar,r))
        x_log.append(np.log10(a_bar(m_bar,r)))

    plt.plot(x_log,y_log)

plt.savefig('figures/all/a_tot_bar_logarithmic')
plt.show()