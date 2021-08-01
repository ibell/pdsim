from PDSim.scroll.scroll_geo import set_scroll_geo
from PDSim.scroll.scroll_geo import geoVals
import PDSim.scroll.symm_scroll_geo as symm

import numpy as np

import pytest

@pytest.mark.parametrize('theta', np.linspace(0, 2*np.pi, 100))
@pytest.mark.parametrize('merge_DDD', [True, False])
@pytest.mark.parametrize('disc', [('2Arc',0.0), ('2Arc', 0.0001), ('ArcLineArc', 0.001)])
def test_fxp_fyp_MOp_sums(theta, merge_DDD, disc, print_values=False):
    """ 
    All the active chambers conntributing to the forces 
    on the orbiting scroll should sum up to zero because of
    a hydrostatic argument.  If the pressure is the same on 
    all sides of the scroll, if the forces and moments
    do not sum to zero, the scroll would spontaneously fly away
    """

    disc_type, disc_r2 = disc

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


@pytest.mark.parametrize('theta', np.linspace(0, 2*np.pi, 11))
@pytest.mark.parametrize('key', ['fx_p','fy_p','M_O_p'])
def test_DD_forces_always_same(theta, key, print_values=True):
    """ 
    All the discharge geometry options must yield the same
    DD force components by a hydrostatic argument
    """
    
    o = []
    disc_options = [('2Arc',0.0), ('2Arc', 0.0001), ('ArcLineArc', 0.001)]
    for i, (disc_type, disc_r2) in enumerate(disc_options):
        geo = geoVals()
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
        o.append(symm.DD_forces(theta, geo, poly=True)[key])
    print(o)
    o = np.array(o)
    is_ok = all(np.abs(o-np.mean(o)) < 1e-8)
    assert(is_ok)

    # assert(size(set(o))==1)

if __name__ == '__main__':
    # print(check_fxp_fyp_MOp_sums(0.4, False, print_values=True))
    pytest.main()