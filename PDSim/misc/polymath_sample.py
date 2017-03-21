# Disable this file, but keep it around in case we want it again
if __name__=='__main__':
    import numpy as np
    from PDSim.misc.polymath import *
    import matplotlib.pyplot as plt

    # An ellipse
    t = np.linspace(0,np.pi,300)
    x1 = np.cos(t)
    y1 = 1.5*np.sin(t)
    p1 = Polygon(x1, y1)

    for offset in [0.0, 0.1, 0.3, 0.4]:
        x2 = np.array([0.0,1.0,1.0,0.0,0.0]) + offset
        y2 = np.array([0.0,0.0,1.0,1.0,0.0]) + offset
        p2 = Polygon(x2,y2)
        
        AND = PolygonOperator(p1, p2).AND()[0]
        
        AND.fill(plt.gca())
        
    plt.gca().set_aspect(1)
    plt.show()