import numpy as np
import scipy.interpolate as sc
import matplotlib.pyplot as plt

x = np.linspace(0, 10, 20)
y = np.sin(x)
tck = sc.splrep(x, y, k=3)
x2 = np.linspace(0, 10, 200)
y2 = sc.splev(x2, tck)

import PDSim.misc.scipylike as psc
spl = psc.splrep(x, y, k=3)
y3 = psc.splev(x2, spl)

plt.plot(x, y, 'o', x2, y2, x2, y3)
plt.show()

print y
for y in y2:
    print y