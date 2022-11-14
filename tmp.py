import numpy as np
from matplotlib import pyplot as plt
x,y,z=[],[],[]
a=[]
for i in range(100):
    x.append(i)
    y.append(i**2)
    z.append(i**3)
    a.append((i**1.5)+100)
# print(a)
plt.figure()
plt.plot(y,z, linestyle = '--')
plt.plot(x,a, color = 'black')
plt.show()