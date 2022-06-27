from scipy import constants
from astropy.units import kg
from astropy.constants import M_sun
from matplotlib import pyplot as plt
import numpy as np
plt.figure()
plt.grid()
# x=[1,2,3,4,5,]
y=[1,4,9,16,25]
x=np.arange(0,5)
plt.plot(x,y)
plt.show()
# plt.scatter(x, y, s=0.5, color='black')
plt.show()