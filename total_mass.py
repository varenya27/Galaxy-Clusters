import csv
from scipy import constants
import numpy as np
from astropy.constants import M_sun
from astropy.units import kg
def M(Th,beta,r_c):
    k = constants.Boltzmann
    mu=0.61
    mp = constants.proton_mass
    G = constants.G
    M_r=((3*beta*k*Th*r_c)/(G*mu*mp))*( 1 / 2)
    # M_r=((3*beta*k*Th)/(G*mu*mp))*((r**2 - r_c**2*ln(r**2+r_c**2))/2 )
    # M_r=((3*beta*k*Th)/(G*mu*mp))*(( r_c**2 - r_c**2*np.log(2*(r_c**2))) + r_c**2*np.log(r_c**2))/2
    return M_r

filename = 'table-1-without_errors.csv'

key,val = [],[]

with open(filename, 'r') as f:
    read = csv.reader(f)
    key = next(read)
    for row in read:
        val.append(row)

# print(key)
# print("------")
# print(val)
names=[]
M_tot=[]
M_tot_solar=[]
for cluster in val:
    names.append(cluster[0])
    Tm=float(cluster[2]) * 11.6*(10**6) #kEv to K
    beta = float(cluster[3])
    r_c = float(cluster[4]) * (3.086*(10**19)) #kpc to m
    M_cluster = M(Tm,beta,r_c)
    M_tot.append(M_cluster)
    M_tot_solar.append("{:e}".format(M_cluster*kg/(M_sun)))

# print(M_tot_solar)
with open('masses.csv','a') as f:
    writer = csv.writer(f,lineterminator='\n')
    for name,mass in zip(names,M_tot_solar):
        writer.writerow([name,mass])