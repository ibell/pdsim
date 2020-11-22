from PDSim.scroll.scroll_geo import set_scroll_geo
from PDSim.scroll.scroll_geo import geoVals
import PDSim.scroll.symm_scroll_geo as symm

import numpy as np

def check_fxp_fyp_MOp_sums(theta, merge_DDD, disc_type, disc_r2, print_values=False):
    """ 
    All the active chambers conntributing to the forces 
    on the orbiting scroll should sum up to zero because of
    a hydrostatic argument.  If the pressure is the same on 
    all sides of the scroll, if the forces and moments
    do not sum to zero, the scroll would spontaneously fly away
    """

    geo = geoVals()
    # set_scroll_geo(100e-6, 2.5, 0.003, 0.005, geo=geo)
    # setDiscGeo(geo, Type='ArcLineArc', r2=0.0008)

    Vdisp = 0.0001048 #[m^3/rev]
    Vratio = 2.2 #[-] 
    t = 0.004 #[m]
    ro = 0.005 #[m]
    phi_i0 = 0.0 #[rad]
    phi_is = 3.141 #[rad]
    phi_os = 0.3 #[rad]    
    #  Set the scroll wrap geometry
    set_scroll_geo(Vdisp, # Vdisp [m^3/rev]
                   Vratio, # Vratio [-]
                   t, # Thickness [m]
                   ro, # Orbiting radius [m]
                   phi_i0 = phi_i0, # [rad]
                   phi_os = phi_os, # [rad]
                   phi_is = phi_is, # [rad]
                   geo=geo)
    symm.setDiscGeo(geo, Type=disc_type, r2=disc_r2)

    Nc = symm.getNc(theta, geo)
    fxp, fyp, MOp = 0, 0, 0
    disc_chambers = ['DDD'] if merge_DDD else ['D1','D2','DD']
    for chamber in ['SA','S1','S2','C1.1','C1.2','C2.1','C2.2'] + disc_chambers:
        if chamber[0:2] in ['C1','C2']:
            func = getattr(symm, chamber[0:2] + "_forces")
            alpha = int(chamber.split('.')[1])
            if alpha > Nc:
                continue
            o = func(theta, alpha, geo)
        else:
            func = getattr(symm, chamber+"_forces")
            o = func(theta, geo)
        if print_values:
            print(chamber, o)
        fxp += o['fx_p']
        fyp += o['fy_p']
        MOp += o['M_O_p']

    if not print_values:
        print(fxp, fyp, MOp)
        assert(abs(fxp) < 1e-16)
        assert(abs(fyp) < 1e-16)
        assert(abs(MOp) < 1e-16)

def test_fxp_fyp_MOpsums():
    for merge_DDD in [True, False]:
        for theta in np.linspace(0, 2*np.pi, 100):
            for disc_type, disc_r2 in [('2Arc',0.0), ('2Arc', 0.0001), ('ArcLineArc', 0.001)]:
                yield check_fxp_fyp_MOp_sums, theta, merge_DDD, disc_type, disc_r2

if __name__ == '__main__':
    # print(check_fxp_fyp_MOp_sums(0.4, False, print_values=True))
    import nose; nose.runmodule()