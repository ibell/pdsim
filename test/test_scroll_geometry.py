from PDSim.scroll.scroll_geo import set_scroll_geo
from PDSim.scroll.symm_scroll_geo import SA_forces
from PDSim.scroll.symm_scroll_geo import setDiscGeo
from PDSim.scroll.scroll_geo import geoVals
import PDSim.scroll.symm_scroll_geo as symm

import numpy as np

def check_fxp_fyp_MOp_sums(theta):
    """ 
    All the active chambers conntributing to the forces 
    on the orbiting scroll should sum up to zero because of
    a hydrostatic argument.  If the pressure is the same on 
    all sides of the scroll, if the forces and and moments
    do not sum to zero, the scroll would spontaneously fly away
    """

    geo = geoVals()
    set_scroll_geo(100e-6, 2.5, 0.003, 0.005, geo=geo)
    setDiscGeo(geo, Type='ArcLineArc', r2=0.0008)

    Nc = symm.getNc(theta, geo)
    fxp, fyp, MOp = 0, 0, 0
    for chamber in ['SA','S1','S2','C1.1','C1.2','C2.1','C2.2','D1','D2','DD']:
        if chamber[0:2] in ['C1','C2']:
            func = getattr(symm, chamber[0:2] + "_forces")
            alpha = int(chamber.split('.')[1])
            if alpha > Nc:
                continue
            o = func(theta, alpha, geo)
        else:
            func = getattr(symm, chamber+"_forces")
            o = func(theta, geo)
        fxp += o['fx_p']
        fyp += o['fy_p']
        MOp += o['M_O_p']
    assert(abs(fxp) < 1e-16)
    assert(abs(fyp) < 1e-16)
    assert(abs(MOp) < 1e-16)

def test_fxp_fyp_MOpsums():
    for theta in np.linspace(0, 2*np.pi, 100):
        yield check_fxp_fyp_MOp_sums, theta

if __name__ == '__main__':

    import nose
    nose.runmodule()