import csv
import numpy as np
from scipy import constants
from scipy.integrate import quad
# n_gas(r) = n_gas(0)*((1+(r/r_c)**2))**(-3*beta/2)
# m_gas = \integral_{0}^{r_c}(n_gas(r)* r^2 * 4\pi)dr
#m_gas =  n_gas(0)


val=[]
with open('table-1-without_errors.csv','r') as f:
    read = csv.reader(f)
    key = next(read)
    for row in read:
        val.append(row)
M_gas =[]
M_error = []
for cluster in val:
    beta = float(cluster[3])
    r_c = float(cluster[4])*(3.08*(10**19))
    n_c = float(cluster[5])*(10**4)
    def M(r):
        n_r=n_c*((1+(r/r_c)**2))**(-3*beta/2)*constants.atomic_mass*0.6
        return n_r*4*np.pi*r**2
    M_cur,err = quad(M,0,r_c)
    M_gas.append(M_cur)
    M_error.append(err)
print(M_gas)
print(M_error)