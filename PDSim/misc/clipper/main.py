'''
Created on 2 oct. 2012 by maxime - Updated by Ian Bell

@author: ibell
Updated 2016/10/20 
@author: Davide Ziviani
'''
from __future__ import division
import sys
from PDSim.misc.clipper import pyclipper
import matplotlib.pyplot as plt
import numpy as np
 
clip = pyclipper.Pyclipper()

square = [[0,0], [1, 0], [1, 1], [0, 1], [0,0]]

x,y = zip(*square)

scale_factor = 10000000

#Rescale the dimensions to go from float to long
scaled_x = [_*scale_factor for _ in x]
scaled_y = [_*scale_factor for _ in y]
clip.subject_polygon([pair for pair in zip(scaled_x,scaled_y)])
plt.plot(x,y)

t = np.linspace(0,2*np.pi,1000)
x = 0.5+0.55*np.cos(t)
y = 0.5+0.55*np.sin(t)
scaled_x = [_*scale_factor for _ in x]
scaled_y = [_*scale_factor for _ in y]
clip.clip_polygon([pair for pair in zip(scaled_x,scaled_y)])
plt.plot(x,y)

sol = clip.execute(pyclipper.INTERSECTION)
for loop in sol:
    scaled_x,scaled_y = zip(*loop)
    x = [_/scale_factor for _ in scaled_x]
    y = [_/scale_factor for _ in scaled_y]
    plt.fill(x,y)

plt.show()