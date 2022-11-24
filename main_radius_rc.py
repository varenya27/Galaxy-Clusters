from scipy import constants
import re
from masses_witherrors import M_total,M_gas,M_stellar,M_stellar2
import numpy as np
from matplotlib import pyplot as plt
import csv
from scipy.optimize import curve_fit
import matplotlib.ticker as tck
from scipy import stats

def line(x, m):
    c=0
    return m*x

k = constants.Boltzmann
mu=0.61
mp = constants.proton_mass
G = constants.G

clusters=[]
with open('tables/table.csv','r') as f:
    read = csv.reader(f)
    key = next(read)
    for row in read:
        clusters.append(row)

for cluster in clusters:
    name=cluster[0]
    print('cluster: '+name)

    x=float((cluster[2][:cluster[2].index("^")])) 
    err_x = (float(re.findall('\{(.*?)\}',cluster[2])[0])-float(re.findall('\{(.*?)\}',cluster[2])[1]))/2
    T=[x*11.6e6,err_x*11.6e6]
    
    x=float((cluster[3][:cluster[3].index("^")])) 
    err_x = (float(re.findall('\{(.*?)\}',cluster[3])[0])-float(re.findall('\{(.*?)\}',cluster[3])[1]))/2
    beta=[x,err_x]
    
    x=float((cluster[4][:cluster[4].index("^")]))
    err_x = (float(re.findall('\{(.*?)\}',cluster[4])[0])-float(re.findall('\{(.*?)\}',cluster[4])[1]))/2
    r_c=[x,err_x]    

    x=float((cluster[5][:cluster[5].index("^")]))
    err_x = (float(re.findall('\{(.*?)\}',cluster[5])[0])-float(re.findall('\{(.*?)\}',cluster[5])[1]))/2
    n_c=[x*(1e4),err_x*(1e4)]   

    x=float((cluster[-2][:cluster[-2].index("^")]))
    err_x = (float(re.findall('\{(.*?)\}',cluster[-2])[0])-float(re.findall('\{(.*?)\}',cluster[-2])[1]))/2
    r500=[x*(3.086e22),err_x*(3.086e22)]   

    x=float((cluster[-1][:cluster[-1].index("^")]))
    err_x = (float(re.findall('\{(.*?)\}',cluster[-1])[0])-float(re.findall('\{(.*?)\}',cluster[-1])[1]))/2
    Th=[x*11.6e6,err_x*11.6e6]    

    t_cool=float((cluster[6][:cluster[6].index("^")]))
    if t_cool<1.3:
        continue
    # print(name)
    # print(name,T,beta,r_c,n_c)
    # threshold=0.2
    # if(T[1]/T[0]>threshold) : continue
    # if(beta[1]/beta[0]>threshold) : continue
    # if(r_c[1]/r_c[0]>threshold) : continue
    # if(n_c[1]/n_c[0]>threshold) : continue

    x,R, err_R=([] for _ in range(3))
    err_tmp,tmp=[],[]
    r_=np.arange(0,r_c[0])
    # for i in range(1,int(r_c[0])):
    for r in r_:
        m_tot,err_m_tot = M_total(r,Th,beta,r_c)
        m_gas, err_m_gas = M_gas(r_c,n_c,beta,r)
        m_stellar, err_m_stellar = M_stellar([m_gas,err_m_gas])
        # m_stellar, err_m_stellar = M_stellar2(M_gas(r_c,n_c,beta,r_c[0]),r500,r)
        
        m_bar,err_m_bar = m_gas+m_stellar, np.sqrt(err_m_gas**2+err_m_stellar**2)
        m_dark,err_m_dark = m_tot-m_bar, np.sqrt(err_m_tot**2+err_m_bar**2)
        if(m_tot*m_gas*m_stellar==0):continue
        if(m_tot/m_bar<=1): 
            start=r
            continue
        # print(err_m_tot,m_tot, err_m_gas,m_gas, err_m_stellar,m_stellar)
        # if(m_stellar<err_m_stellar): print('starf')
        # if(m_gas<err_m_gas): print('gas f')
        # if(m_tot<err_m_tot): print('tot f')
        # if(m_bar<err_m_bar): print('bar f')
        # if(m_dark<err_m_dark): 
        #     print('darkf {:e} {:e} {:e} {:e}'.format(m_tot,m_gas, m_stellar,err_m_tot))

        x.append(r)
        mratio=m_dark/m_bar
        # print("{:e} {:e} {:e} {:e}".format(m_dark,err_m_dark,m_stellar,err_m_stellar))
        # print(err_m_bar/m_bar, err_m_gas/m_gas,err_m_stellar/m_stellar,)
        err_mratio = np.sqrt( (err_m_dark/m_bar)**2+(m_dark*err_m_bar/(m_bar**2))**2 )
        R.append(mratio)
        err_R.append(err_mratio)
        tmp.append(m_stellar)
        err_tmp.append(err_m_stellar)
        print('{:e} {:e} {:e}'.format(m_tot,m_gas, m_stellar,))
    print('-----')
    N=len(R)
    N=int(N/10)+1

    param, param_cov = curve_fit(line, x, R, sigma=err_R, absolute_sigma=True)
    perr = np.sqrt(np.diag(param_cov))
    R_result = param[0]*np.array(x)#+param[1]
    # print(R,err_R)
    delta = np.array(R[::N]) - R_result[::N]
    chi_sq = np.sum(delta**2/np.array(err_R[::N])**2)
    df = len(R_result[::N])
    p=  1 - stats.chi2.cdf(chi_sq,df)
    # print(p)
    fig, ax = plt.subplots()
    ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.patch.set_edgecolor('black')  
    ax.patch.set_linewidth('1')
    plt.title('Cluster '+name+' | m='+str(round(param[0],2))+' | $\chi^2$= '+str(round(chi_sq,2)))
    plt.xlabel('r (kpc)')
    plt.ylabel('$M_{dark}/M_{bar}$')
    plt.errorbar(x[::N], R[::N],err_R[::N],linestyle='',fmt='h',color='black' ,alpha=0.5, capsize=2.0, zorder=0,label ="Data")
    # plt.errorbar(x[::N], tmp[::N],err_tmp[::N],linestyle='',fmt='h',color='blue' ,alpha=0.5, capsize=2.0, zorder=0,label ="$M_{star}$")
    plt.plot(x, R_result, '--', label ="Best Fit",color='orange')
    plt.legend()
    # plt.savefig('figures/'+name+'.png')
    plt.show()
    plt.close()


    with open('chi.txt','a') as f:
        x = name+' : '+str(round(chi_sq,2))+' : '+str(p)+'\n'
        f.write(x)
