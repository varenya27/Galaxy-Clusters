import re
import numpy as np
import csv
from scipy import constants
from scipy.optimize import curve_fit
from scipy import stats
import matplotlib.ticker as tck
from matplotlib import pyplot as plt
from matplotlib import gridspec
from masses_witherrors import M_total,M_gas,M_stellar,M_stellar2

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
count=0
for cluster in clusters:
    name=cluster[0]
    # print('cluster: '+name)

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
    r500=[x,err_x]   

    x=float((cluster[-1][:cluster[-1].index("^")]))
    err_x = (float(re.findall('\{(.*?)\}',cluster[-1])[0])-float(re.findall('\{(.*?)\}',cluster[-1])[1]))/2
    Th=[x*11.6e6,err_x*11.6e6]    

    t_cool=float((cluster[6][:cluster[6].index("^")]))
    if t_cool<1.3:
        continue
    count+=1
    # print(name)
    # print(name,T,beta,r_c,n_c)
    # threshold=0.2
    # if(T[1]/T[0]>threshold) : continue
    # if(beta[1]/beta[0]>threshold) : continue
    # if(r_c[1]/r_c[0]>threshold) : continue
    # if(n_c[1]/n_c[0]>threshold) : continue

    x,R, err_R,Mdark,Mbar,Mtot=([] for _ in range(6))
    err_tmp,tmp=[],[]
    rc_list=np.arange(2,r_c[0])
    # rang=range(1,int(r500[0]))#[::int(r500[0]/20)]
    # for r in range(1,int(r500[0])):
    r500_list=np.arange(2,r500[0]*1000)
    # print(name,r500_list)
    for r in r500_list:
    # for r in rc_list:
        # r/=1000
        # print(r)
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
        Mdark.append(m_dark)
        Mtot.append(m_tot)
        Mbar.append(m_bar)
    # print('{:f} {:f} {:f} {:f}'.format(m_tot/1e14,err_m_tot/1e14,m_gas/1e13, err_m_gas/1e13))
    # print(name,' -----')
    N=len(R)
    N=int(N/10)+1
    # quit()
    param, param_cov = curve_fit(line, x, R, sigma=err_R, absolute_sigma=True)
    perr = np.sqrt(np.diag(param_cov))
    R_result = param[0]*np.array(x)#+param[1]
    # print(R,err_R)
    delta = np.array(R[::N]) - R_result[::N]
    chi_sq = np.sum(delta**2/np.array(err_R[::N])**2)
    df = len(R_result[::N])
    p=  1 - stats.chi2.cdf(chi_sq,df)
    # print(p)

    fig=plt.figure()
    fig.set_figheight(7)
    fig.set_figwidth(6)
    spec = gridspec.GridSpec(ncols=1, nrows=2,hspace=0.5, height_ratios=[2, 1])
    ax1= fig.add_subplot(spec[0])
    ax1.set_title('Cluster '+name+' | $\chi^2$= '+str(round(chi_sq,2)))
    ax1.set_xlabel('r (Mpc)')
    ax1.set_ylabel('$M_{dark}/M_{bar}$')
    ax1.yaxis.set_minor_locator(tck.AutoMinorLocator())
    ax1.xaxis.set_minor_locator(tck.AutoMinorLocator(10))
    ax1.patch.set_edgecolor('black')  
    ax1.patch.set_linewidth('1')
    ax1.tick_params(axis="x", which="minor", direction="in",top=True,bottom=True,)
    ax1.tick_params(axis="y", which="minor", direction="in",left=True,right=True,)
    ax1.tick_params(axis="x", which="major", direction="inout",top=True,bottom=True,)
    ax1.tick_params(axis="y", which="major", direction="inout",left=True,right=True,)

    ax1.errorbar(x[::50], R[::50],err_R[::50],linestyle='',fmt='h',color='black' ,alpha=0.5, capsize=2.0, zorder=0,label ="Data")
    # plt.errorbar(x[::25], R[::25],linestyle='',fmt='o',color='black' ,alpha=0.5, capsize=2.0, zorder=0,)
    # plt.errorbar(x[::N], tmp[::N],err_tmp[::N],linestyle='',fmt='h',color='blue' ,alpha=0.5, capsize=2.0, zorder=0,label ="$M_{star}$")
    ax1.plot(x, R_result, '--', label ="Best Fit",color='darkred')
    ax1.legend(loc='lower right')
    # ax1.legend(loc='lower right')
    # plt.figure()
    
    ax2 = fig.add_subplot(spec[1])
    ax2.plot(x, Mdark, '-',linewidth=3, label ="Dark Mass",)
    ax2.plot(x, Mbar, '-',linewidth=3, label ="Baryonic Mass",)
    ax2.plot(x, Mtot, '-',linewidth=3, label ="Total Mass",)    
    ax2.set_title('Mass Profiles')
    ax2.set_xlabel('r (Mpc)')
    ax2.set_ylabel('$M (M_{\odot})$')
    ax2.yaxis.set_minor_locator(tck.AutoMinorLocator())
    ax2.xaxis.set_minor_locator(tck.AutoMinorLocator(10))
    ax2.patch.set_edgecolor('black')  
    ax2.patch.set_linewidth('1')    
    ax2.tick_params(axis="x", which="minor", direction="in",top=True,bottom=True,)
    ax2.tick_params(axis="y", which="minor", direction="in",left=True,right=True,)
    ax2.tick_params(axis="x", which="major", direction="inout",top=True,bottom=True,)
    ax2.tick_params(axis="y", which="major", direction="inout",left=True,right=True,)

    ax2.legend(loc='best')


    # fig.tight_layout()
    plt.savefig('figures/'+name+'.png')
    # plt.show()
    plt.close()


    with open('chi.txt','a') as f:
        x = name+' : '+str(round(chi_sq,2))+' : '+str(p)+'\n'
        f.write(x)
print('total number of clusters:',count)