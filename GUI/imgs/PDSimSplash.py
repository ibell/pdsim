import matplotlib.pyplot as plt
import numpy as np
from math import pi
import PDSim
plt.figure(figsize=(6,3))

w=1
h=1
x=np.linspace(-w,w,1000)
y=np.linspace(-h,h,1000)
X,Y=np.meshgrid(x,y)
Z=np.cos(X*pi/(2*w))**(0.4)*np.cos(Y*pi/(2*h))**(0.4)
plt.gray()
plt.contourf(X,Y,Z,100,alpha=1.0)
plt.text(0,h/2.0,'RecipGUI',ha='center',va='center',size=40,color='black',alpha=0.6)
plt.text(0+0.01,h/2.0+0.02,'RecipGUI',ha='center',va='center',size=40,color='red')

version=PDSim.__version__
plt.text(0,-h/2.0,'Powered by PDSim v. '+version,size=20,ha='center',va='center')
plt.text(0,-h/1.25,'by Ian Bell',size=12,ha='center',va='center')
plt.gca().axis('off')
plt.savefig('RecipGUI.png',dpi=300)
plt.show()