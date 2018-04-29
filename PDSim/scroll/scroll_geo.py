# -*- coding: utf-8 -*-
from __future__ import absolute_import

from .common_scroll_geo import *
from .symm_scroll_geo import *

from math import pi

# If scipy is available, use its interpolation and optimization functions, otherwise, 
# use our implementation (for packaging purposes mostly)
try:
    from scipy.optimize import fsolve
except ImportError:
    from PDSim.misc.solvers import MultiDimNewtRaph as fsolve

def set_scroll_geo(Vdisp,Vratio,Thickness,OrbitingRadius,phi_i0=0.0,phi_os=0.3, phi_is = pi, geo = None):
    """
    Provide the following parameters.  The rest will be calculated by the geometry code
    
    ==============  ===================================================================================
    Vdisp           Displacement in compressor mode [m^3]
    Vratio          Volume ratio (compression chambers at discharge angle / displacement volume) [-]
    Thickness       Thickness of scroll wrap [m]
    OrbitingRadius  Orbiting radius of the orbiting scroll [m]
    ==============  ===================================================================================
    
    Optional parameters are 
    
    phi_i0
    phi_os
    phi_is
    """
    
    ## Determine the geometry by using the imposed parameters for the scroll wraps
    def f(x,phi_i0,phi_os,Vdisp_goal,Vratio_goal,t_goal,ro_goal):
        phi_ie=x[0]
        phi_o0=x[1]
        hs=x[2]
        rb=x[3]
        t=rb*(phi_i0-phi_o0)
        ro=rb*pi-t
        Vdisp=-2*pi*hs*rb*ro*(3*pi-2*phi_ie+phi_i0+phi_o0)
        Vratio=(3*pi-2*phi_ie+phi_i0+phi_o0)/(-2*phi_os-3*pi+phi_i0+phi_o0)

        r1=Vdisp-Vdisp_goal
        r2=Vratio-Vratio_goal
        r3=t-t_goal
        r4=ro-ro_goal
        return [r1,r2,r3,r4]
    
    phi_ie,phi_o0,hs,rb = fsolve(f,[20,1.3,0.03,0.003],args=(phi_i0,phi_os,Vdisp,Vratio,Thickness,OrbitingRadius))
    phi_oe=phi_ie

    # Return the values
    if geo is None:
        geo = geoVals()
    geo.h=hs
    geo.rb=rb
    geo.phi_oi0=geo.phi_fi0=phi_i0
    geo.phi_ois=geo.phi_fis=phi_is
    geo.phi_oie=geo.phi_fie=phi_ie
    geo.phi_oo0=geo.phi_fo0=phi_o0
    geo.phi_oos=geo.phi_fos=phi_os
    geo.phi_ooe=geo.phi_foe=phi_oe
    geo.ro=rb*pi-Thickness
    geo.t=Thickness
    return geo