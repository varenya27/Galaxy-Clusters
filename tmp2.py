from scipy import constants
import re
from masses_witherrors import M_total,M_gas,M_stellar,M_stellar2
import numpy as np
from matplotlib import pyplot as plt
import csv
from scipy.optimize import curve_fit
import matplotlib.ticker as tck
from scipy import stats
x=np.linspace(0,10,1000)
R=x*np.exp(-x)
err_R = R/10
R_result = 0.01*x

Mdark=np.exp(-2*x)
Mbar=np.exp(-3*x)
Mtot=np.exp(-1*x)
name='test'
chi_sq=5.23
if(1):

    # fig, ax = plt.subplots(nrows=2, ncols=1)
    fig=plt.figure()
    # set height of each subplot as 8
    fig.set_figheight(7)
    
    # set width of each subplot as 8
    fig.set_figwidth(6)
    spec = gridspec.GridSpec(ncols=1, nrows=2,
                         hspace=0.5, height_ratios=[2, 1])
    # fig.set_title('asdf')
    ax1= fig.add_subplot(spec[0])
    ax1.set_title('Cluster '+name+' | $\chi^2$= '+str(round(chi_sq,2)))
    ax1.set_xlabel('r (Mpc)')
    ax1.set_ylabel('$M_{dark}/M_{bar}$')
    ax1.yaxis.set_minor_locator(tck.AutoMinorLocator())
    ax1.xaxis.set_minor_locator(tck.AutoMinorLocator(10))
    ax1.patch.set_edgecolor('black')  
    ax1.patch.set_linewidth('1')

    ax1.errorbar(x[::50], R[::50],err_R[::50],linestyle='',fmt='h',color='black' ,alpha=0.5, capsize=2.0, zorder=0,label ="Data")
    # plt.errorbar(x[::25], R[::25],linestyle='',fmt='o',color='black' ,alpha=0.5, capsize=2.0, zorder=0,)
    # plt.errorbar(x[::N], tmp[::N],err_tmp[::N],linestyle='',fmt='h',color='blue' ,alpha=0.5, capsize=2.0, zorder=0,label ="$M_{star}$")
    ax1.plot(x, R_result, '--', label ="Best Fit",color='darkred')
    ax1.legend()
    # plt.figure()
    
    ax2 = fig.add_subplot(spec[1])
    ax2.plot(x, Mdark, '-',linewidth=3, label ="Dark Mass",)
    ax2.plot(x, Mbar, '-',linewidth=3, label ="Baryonic Mass",)
    ax2.plot(x, Mtot, '-',linewidth=3, label ="Total Mass",)    
    ax2.set_xlabel('r (Mpc)')
    ax2.set_ylabel('$M (M_{\odot})$')
    ax2.yaxis.set_minor_locator(tck.AutoMinorLocator())
    ax2.xaxis.set_minor_locator(tck.AutoMinorLocator(10))
    ax2.patch.set_edgecolor('black')  
    ax2.patch.set_linewidth('1')    
    
    ax2.legend()


    # fig.tight_layout()
    # plt.savefig('figures/'+name+'.png')
    plt.show()
    plt.close()

