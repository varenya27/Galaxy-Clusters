import csv
import numpy as np
from scipy import constants
from astropy.constants import M_sun
from astropy.units import kg

def M(m_tot):
    return 4*(10**12)*(((m_tot)/(5.7*(10**13)))**0.6)

M_stellar =[]
val=[]
with open ('masses.csv','r') as f:
    read = csv.reader(f)
    key = next(read)
    for row in read:
        m_tot=float(row[1])
        M_stellar.append("{:e}".format(M(m_tot)))
print(M_stellar)

with open('masses2.csv','r') as i:
    with open('masses3.csv','w') as o:
        writer = csv.writer(o,lineterminator='\n')
        data=[]
        data.append(next(csv.reader(i)))
        for row,star in zip(csv.reader(i),M_stellar):
            row.append(star)
            data.append(row)
        writer.writerows(data)