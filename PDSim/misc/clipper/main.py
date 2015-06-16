'''
Created on 2 oct. 2012 by maxime - Updated by Ian Bell
Update 16 June 2015 by Davide Ziviani

@author: ibell
@contributor: dziviani
'''
from __future__ import division
import sys

from PDSim.misc.clipper import pyclipper
import matplotlib.pyplot as plt
import numpy as np


clip = pyclipper.Pyclipper()
scale_factor = 10000000


def Test(clip,scale_factor):
    """
    Clipping polygon: Square 
    """
    
    square = [[0,0], [1, 0], [1, 1], [0, 1], [0,0]]
    x1,y1 = zip(*square)
    scale_factor = 10000000
    #Rescale the dimensions to go from float to long
    scaled_x = [_*scale_factor for _ in x1]
    scaled_y = [_*scale_factor for _ in y1]
    clip.clip_polygon([pair for pair in zip(scaled_x,scaled_y)])
    #clip.add_polygon([pair for pair in zip(scaled_x,scaled_y)])


    """
    Subject polygon: Ellipse
    """
    scale_factor = 10000000
    t = np.linspace(0,2*np.pi,1000)
    x2 = 0.5+0.55*np.cos(t)
    y2 = 0.5+0.55*np.sin(t)
    scaled_x = [_*scale_factor for _ in x2]
    scaled_y = [_*scale_factor for _ in y2]
    clip.subject_polygon([pair for pair in zip(scaled_x,scaled_y)])
    #clip.sub_polygon([pair for pair in zip(scaled_x,scaled_y)])
    
    return x1,y1,x2,y2,clip



def polygon_intersect(clip,scale_factor):
    """
    Methods:
    INTERSECTION = ctIntersection
    DIFFERENCE = ctDifference
    UNION = ctUnion
    XOR = ctXor
    """
    sol = clip.execute(pyclipper.INTERSECTION)
    for loop in sol:
        scaled_x,scaled_y = zip(*loop)
        x = [_/scale_factor for _ in scaled_x]
        y = [_/scale_factor for _ in scaled_y]
        plt.fill(x,y,'r',alpha=0.3)




x1,y1,x2,y2,clip = Test(clip,scale_factor)

polygon_intersect(clip,scale_factor)


plt.plot(x1,y1,'k-',linewidth = 1.5)    
plt.plot(x2,y2,'b-',linewidth = 1.5)

plt.xlim(-1,2)
plt.ylim(-1,2)
plt.show()