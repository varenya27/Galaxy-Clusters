from scipy import constants
import re
from masses_witherrors import M_total,M_gas,M_stellar
import numpy as np
from matplotlib import pyplot as plt
import csv
from scipy.optimize import curve_fit

def line(x, m):
    c=0
    return m*x

k = constants.Boltzmann
mu=0.61
mp = constants.proton_mass
G = constants.G
n_iter=1

clusters=[]
with open('tables/table-1.csv','r') as f:
    read = csv.reader(f)
    key = next(read)
    for row in read:
        clusters.append(row)

for cluster in clusters:
    name=cluster[0]

    x=float((cluster[2][:cluster[2].index("^")])) 
    err_x = (float(re.findall('\{(.*?)\}',cluster[2])[0])-float(re.findall('\{(.*?)\}',cluster[2])[1]))/2
    T=[x*11.6e6,err_x*11.6e6]
    
    x=float((cluster[3][:cluster[3].index("^")])) 
    err_x = (float(re.findall('\{(.*?)\}',cluster[3])[0])-float(re.findall('\{(.*?)\}',cluster[3])[1]))/2
    beta=[x,err_x]
    
    x=float((cluster[4][:cluster[4].index("^")]))
    err_x = (float(re.findall('\{(.*?)\}',cluster[4])[0])-float(re.findall('\{(.*?)\}',cluster[4])[1]))/2
    r_c=[x,err_x]    

    x=float((cluster[5][:cluster[5].index("^")]))*(1e4)
    err_x = (float(re.findall('\{(.*?)\}',cluster[5])[0])-float(re.findall('\{(.*?)\}',cluster[5])[1]))/2*(1e4)
    n_c=[x,err_x]    

    t_cool=float((cluster[6][:cluster[6].index("^")]))
    if t_cool<1.3:
        continue

    x,R, err_R=([] for _ in range(3))

    for i in range(0,n_iter*int(r_c[0])):
        r=i/n_iter
        m_tot,err_m_tot = M_total(r,T,beta,r_c)
        m_gas, err_m_gas = M_gas(r_c,n_c,beta,r)
        m_stellar, err_m_stellar = M_stellar([m_gas,err_m_gas])
        m_bar,err_m_bar = m_gas+m_stellar, err_m_gas+err_m_stellar
        m_dark,err_m_dark = m_tot-m_bar, err_m_tot-err_m_bar
        if(m_tot*m_gas*m_stellar==0):continue
        if(m_tot/m_bar<=1): 
            start=i
            continue
        # print(err_m_tot,m_tot, err_m_gas,m_gas, err_m_stellar,m_stellar)
        # if(m_stellar>err_m_stellar): print('yes')
        # else: print('no')
        x.append(i)
        mratio=m_dark/m_bar
        # print("{:e} {:e} {:e} {:e}".format(m_dark,err_m_dark,m_gas,err_m_gas))
        err_mratio = err_m_dark/m_bar-m_dark*err_m_bar/(m_bar**2)
        R.append(mratio)
        err_R.append(err_mratio)
        
    # print(R,err_R)


    param, param_cov = curve_fit(line, x[1:], R[1:], sigma=err_R[1:], absolute_sigma=True)
    perr = np.sqrt(np.diag(param_cov))
    R_result = param[0]*np.array(x)#+param[1]
    # print(R,err_R)
    delta = np.array(R) - R_result
    chi_sq = np.sum(delta[1:]**2/np.array(err_R[1:])**2)
    print("parameters: ",param)
    plt.figure()
    plt.grid()
    plt.title('Cluster '+name)
    plt.xlabel('r (kpc)')
    # plt.errorbar(x, R,err_R)
    plt.scatter(x[::5], R[::5],3, label ="data")
    plt.plot(x, R_result, '--', label ="fit",color='orange')
    plt.legend()
    plt.show()
    print(chi_sq)
    break