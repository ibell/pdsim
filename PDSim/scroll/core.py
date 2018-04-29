from __future__ import division, print_function, absolute_import
##--- package imports
from PDSim.misc.datatypes import arraym
from PDSim.core.containers import ControlVolume
from PDSim.flow.flow import FlowPath
from PDSim.core.core import PDSimCore
from PDSim.flow import flow_models

from PDSim.core.bearings import journal_bearing
from PDSim.scroll import scroll_geo, symm_scroll_geo
from ._scroll import _Scroll
from PDSim.scroll.scroll_geo import set_scroll_geo

##--- non-package imports
import warnings

from CoolProp import State
from math import pi
import numpy as np
import copy
import types

# If scipy is available, use its interpolation and optimization functions, otherwise, 
# use our implementation (for packaging purposes mostly)
try:
    import scipy.interpolate as interp
    import scipy.optimize as optimize
    from scipy.optimize import fsolve
except ImportError:
    import PDSim.misc.scipylike as interp
    import PDSim.misc.solvers as optimize
    from PDSim.misc.solvers import MultiDimNewtRaph as fsolve

import matplotlib.pyplot as plt
import subprocess
import glob
import os

class struct(object):
    pass

class Port(object):
    
    #: Involute angle of the involute used to locate this port
    phi = 3.14159
    
    #: The code for the involute used to locate this point: 'i' or 'o'
    involute = 'i'
    
    #: Distance away from the involute
    offset = 0.001
    
    #: Diameter of the port
    D = 0.0005
    
    #: The x coordinate of the center of the port
    x0 = None
    
    #: The y coordinate of the center of the port
    y0 = None
    
    #: The array of crank angles used to calculate the free area of the port
    theta = None
    
    #: A dictionary with keys of each partner CV, and values equal to the 
    #: free area of the port with that control volume
    area_dict = None

class Scroll(PDSimCore, _Scroll):
    """
    This is a python class that implements functionality for a scroll compressor
    
    It is inherited from the PDSimCore class
    
    This class is only to be used for symmetric scroll wraps, for asymmetric 
    scroll wraps, use the AsymmetricScroll class
    """
    
    def __init__(self):
        PDSimCore.__init__(self)
        
        ## Define the geometry structure
        self.geo=scroll_geo.geoVals()

        ## Structure for virtual sensors        
        self.sensors = struct()
        
        ## Set flags
        self.__Setscroll_geo__=False
        self.__SetDiscGeo__=False
        self.__before_discharge__=False #Step bridging theta_d

    def __getstate__(self):
        """
        A function for preparing class instance for pickling
         
        Combine the dictionaries from the _Scroll base class and the Scroll
        class when pickling
        """
        py_dict = self.__dict__.copy()
        py_dict.update(_Scroll.__cdict__(self))
        return py_dict

    def __setstate__(self, d):
        """
        A function for unpacking class instance for unpickling
        """
        for k,v in d.iteritems():
            setattr(self,k,v)
        
    def INTERPOLATING_NOZZLE_FLOW(self, 
                                FlowPath, 
                                X_d = 1.0,
                                X_d_backflow = 0.8,
                                upstream_key = 'EVICIV',
                                A_interpolator = None,
                                DP_floor = 1e-10):
        """
        A generic isentropic nozzle flow model wrapper with the added consideration
        that if the pressure drop is below the floor value, there is no flow.
        This code was originally added to handle the case of the injection line 
        where there is no flow out of the injection which greatly increases the 
        numerical stiffness.
        
        Furthermore, the area is determined through the use of the spline
        interpolator
        
        This function also implements the use of spline interpolation to calculate
        the area between the EVI port and the control volume
        
        Parameters
        ----------
        FlowPath : FlowPath instance
            A fully-instantiated flow path model
        X_d : float
            Flow coefficient when the flow goes from ``upstream_key`` to the downstream key
        X_d_backflow : float
            Flow coefficient when the flow goes from downstream key to  ``upstream_key``
        upstream_key : string
            Key for the side of the flow path that is considered to be "upstream"
        A_interpolator : float
            throat area for isentropic nozzle model [:math:`m^2`]
        DP_floor: float
            The minimum pressure drop [kPa]
            
        Returns
        -------
        mdot : float
            The mass flow through the flow path [kg/s]
        """
        FlowPath.A = interp.splev(self.theta, A_interpolator)
        
        try:
            if FlowPath.State_up.p - FlowPath.State_down.p > DP_floor and abs(FlowPath.A) > 1e-15:
                if FlowPath.key_up == upstream_key:
                    _X_d = X_d #Normal flow into CV from EVI
                else:
                    _X_d = X_d_backflow #backflow from CV to EVI
                mdot = _X_d*flow_models.IsentropicNozzle(FlowPath.A,
                                                         FlowPath.State_up,
                                                         FlowPath.State_down)
                return mdot
            else:
                return 0.0
        except ZeroDivisionError:
            return 0.0
            
    def calculate_port_areas(self):
        """ 
        Calculate the area between a port on the fixed scroll and all of the 
        control volumes.

        This port could be a port for injection, a pressure tap for dynamic
        pressure measurement, etc.
        
        This function iterates over the ports in ``self.fixed_scroll_ports`` and
        loads the variables ``theta`` and ``area_dict`` for each port
        """
        
        # Iterate over the ports
        for port in self.fixed_scroll_ports:
            
            #  Make sure it is an Port instance
            assert (isinstance(port,Port))
              
            #  Get the reference point on the scroll wrap
            if port.involute == 'i':
                #  Point on the scroll involute
                x, y = scroll_geo.coords_inv(port.phi, self.geo, 0, 'fi')
                #  Unit normal vector components
                nx, ny = scroll_geo.coords_norm(port.phi, self.geo, 0, 'fi')
            elif port.involute == 'o':
                #  Point on the scroll involute
                x, y = scroll_geo.coords_inv(port.phi, self.geo, 0, 'fo')
                #  Unit normal vector components
                nx, ny = scroll_geo.coords_norm(port.phi, self.geo, 0, 'fo')
            else:
                raise ValueError('port involute[{0:s}] must be one of "i" or "o"'.format(port.involute))
            
            #  Normal direction points towards the scroll wrap, take the opposite 
            #  direction to locate the center of the port
            port.x0 = x - port.offset*nx
            port.y0 = y - port.offset*ny
            
            #  The coordinates for the center of the port
            t = np.linspace(0, 2*pi)
            xport = port.x0 + port.D/2.0*np.cos(t)
            yport = port.y0 + port.D/2.0*np.sin(t)
            
            #  Actually use the clipper library to carry out the intersection
            #  of the port with all of the control volumes
            theta_area, area_dict = self.poly_intersection_with_cvs(xport, yport, 100)
            
            #  Save the values
            port.area_dict = area_dict
            port.theta = theta_area
            
#            #  Plot them
#            for k, A in area_dict.iteritems():
#                plt.plot(theta_area, A*1e6, label = k)
#            
#            plt.legend()
#            plt.xlabel('Crank angle [rad]')
#            plt.ylabel('Area [mm$^2$]')
#            plt.savefig('Aport.png')
#            plt.show()  
            
    def get_discharge_port_blockage_poly(self, theta):
        """
        Get all the polygons associated with the control volumes that could in principle
        be connected with the discharge port
        """
        xdd, ydd = scroll_geo.CVcoords('dd',self.geo,theta)
        xd1, yd1 = scroll_geo.CVcoords('d1',self.geo,theta)
        Ncmax = scroll_geo.nC_Max(self.geo)
        Nc = scroll_geo.getNc(theta, self.geo)
        
        if Nc == Ncmax and Nc > 0:
            # Inner-most compression chamber in use
            xc1_N, yc1_N = scroll_geo.CVcoords('c1.+'+str(Ncmax), self.geo, theta)
            xc2_N, yc2_N = scroll_geo.CVcoords('c2.+'+str(Ncmax), self.geo, theta)
        elif Nc == Ncmax and Nc == 0:
            # Suction chamber 
            xc1_N, yc1_N = scroll_geo.CVcoords('s1', self.geo, theta)
            xc2_N, yc2_N = scroll_geo.CVcoords('s2', self.geo, theta)
        else:
            xc1_N, yc1_N = None, None
            xc2_N, yc2_N = None, None
            
        if Nc == Ncmax-1 and Ncmax > 0:
            xc1_Nm1, yc1_Nm1 = scroll_geo.CVcoords('c1.+'+str(Ncmax-1), self.geo, theta)
        else:
            xc1_Nm1, yc1_Nm1 = None, None
            
        return dict(xdd = xdd, ydd = ydd,
                    xd1 = xd1, yd1 = yd1,
                    xc1_N = xc1_N, yc1_N = yc1_N,
                    xc1_Nm1 = xc1_Nm1, yc1_Nm1 = yc1_Nm1
                    )
                
    def cache_discharge_port_blockage(self, xport = None, yport = None, plot = False, N = 100):
        """
        Precalculate the discharge port blockage using the clipper polygon math module
        
        This is computationally quite expensive, which is why it is done only once
        and then interpolated within.
        
        The port position defaults to be equal to the center of arc1 in the discharge 
        region with 90% the radius of arc1
        
        Parameters
        ----------
        xport : list or 1D numpy array of x coordinates, optional
            The x coordinates for the port
        yport : list or 1D numpy array of y coordinates, optional
            The y coordinates for the port
        plot  : bool, optional
            Whether or not to generate plots for each crank angle
        N : int, optional
            How many steps to include over one rotation
        """
        
        if plot:
            print('plotting of disc port blockage is on')
            
        from PDSim.misc.clipper import pyclipper
        
        scale_factor = 1000000000
        
        if xport is None and yport is None:
            warnings.warn('xport and yport not provided, defaulting back to circular discharge port; should be stored in self.geo.xvec_disc_port and self.geo.yvec_disc_port')
            xport = self.geo.xa_arc1 + 0.9*self.geo.ra_arc1*np.cos(np.linspace(0,2*pi,100))
            yport = self.geo.ya_arc1 + 0.9*self.geo.ra_arc1*np.sin(np.linspace(0,2*pi,100))
        
        print('caching discharge port blockage, please wait...', end=' ')
        # Scale the floating points to long integers
        scaled_xport = xport*scale_factor
        scaled_yport = yport*scale_factor
        
        xinvi, yinvi = scroll_geo.coords_inv(np.linspace(self.geo.phi_fis+2*np.pi,self.geo.phi_fis,75), self.geo, 0, flag="fi")
        xarc1 = self.geo.xa_arc1 + self.geo.ra_arc1*np.cos(np.linspace(self.geo.t2_arc1,self.geo.t1_arc1,75))
        yarc1 = self.geo.ya_arc1 + self.geo.ra_arc1*np.sin(np.linspace(self.geo.t2_arc1,self.geo.t1_arc1,75))
        xarc2 = self.geo.xa_arc2 + self.geo.ra_arc2*np.cos(np.linspace(self.geo.t1_arc2,self.geo.t2_arc2,75))
        yarc2 = self.geo.ya_arc2 + self.geo.ra_arc2*np.sin(np.linspace(self.geo.t1_arc2,self.geo.t2_arc2,75))
        xinvo, yinvo = scroll_geo.coords_inv(np.linspace(self.geo.phi_fos,self.geo.phi_fis+2*np.pi,75), self.geo, 0, flag="fo")
        
        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        t, A, Add, Ad1, Ac1_N, Ac1_Nm1 = [], [], [], [], [], []
        for i,theta in enumerate(np.linspace(0, 2*pi, 100)):
            
            THETA = self.geo.phi_fie-pi/2.0-theta

            # Get the polygons from the polygon generation function            
            poly = self.get_discharge_port_blockage_poly(theta)
            xdd, ydd = poly['xdd'],poly['ydd']
            xd1, yd1 = poly['xd1'],poly['yd1']
            xc1_N, yc1_N = poly['xc1_N'],poly['yc1_N']
            xc1_Nm1, yc1_Nm1 = poly['xc1_Nm1'],poly['yc1_Nm1']
            
            scaled_xdd = xdd*scale_factor
            scaled_ydd = ydd*scale_factor
            scaled_xd1 = xd1*scale_factor
            scaled_yd1 = yd1*scale_factor
            if xc1_N is not None:
                scaled_xc1_N = xc1_N*scale_factor
                scaled_yc1_N = yc1_N*scale_factor
            else:
                scaled_xc1_N = None
                scaled_yc1_N = None
            if xc1_Nm1 is not None:
                scaled_xc1_Nm1 = xc1_Nm1*scale_factor
                scaled_yc1_Nm1 = yc1_Nm1*scale_factor
            else:
                scaled_xc1_Nm1 = None
                scaled_yc1_Nm1 = None
            
            if plot:
                ax.cla()
        
            xscroll = -np.r_[xinvi,xarc1,xarc2,xinvo,xinvi[0]]+self.geo.ro*np.cos(THETA)
            yscroll = -np.r_[yinvi,yarc1,yarc2,yinvo,yinvi[0]]+self.geo.ro*np.sin(THETA)
            
            xscroll = xscroll[::-1]
            yscroll = yscroll[::-1]
            
            scaled_xscroll = xscroll*scale_factor
            scaled_yscroll = yscroll*scale_factor
            
            clip = pyclipper.Pyclipper()
            clip.subject_polygon([pair for pair in zip(scaled_xport, scaled_yport)])
            clip.clip_polygon([pair for pair in zip(scaled_xscroll, scaled_yscroll)])
            soln = clip.execute(pyclipper.DIFFERENCE)
            
            if plot:
                ax.plot(scaled_xport, scaled_yport)
                ax.plot(scaled_xscroll, scaled_yscroll)
                ax.fill(scaled_xdd, scaled_ydd, alpha = 0.1, zorder = 1000)
                ax.fill(scaled_xd1, scaled_yd1, alpha = 0.1, zorder = 1000)
                if xc1_N is not None:
                    ax.fill(scaled_xc1_N, scaled_yc1_N, alpha = 0.1, zorder = 1000)
                if xc1_Nm1 is not None:
                    ax.fill(scaled_xc1_Nm1, scaled_yc1_Nm1, alpha = 0.1, zorder = 1000)
                
            for loop in soln:
                scaled_x, scaled_y = zip(*loop)
                if plot:
                    ax.fill(scaled_x, scaled_y)
            
            #  The total flow area unblocked by the scroll wrap    
            Atotal = sum([pyclipper.area(loop) for loop in soln])/scale_factor**2
            
            def calculate_area(scaled_xcv,scaled_ycv):
                if scaled_xcv is None:
                    return 0.0
                    
                A_CV = 0
                for loop in soln:
                    scaled_x, scaled_y = zip(*loop)
                    if plot:
                        ax.fill(scaled_x, scaled_y, 'r')
                    
                    #  Now see if it has overlap with the DD chamber    
                    clip = pyclipper.Pyclipper()
                    clip.subject_polygon([pair for pair in zip(scaled_x, scaled_y)])
                    clip.clip_polygon([pair for pair in zip(scaled_xcv, scaled_ycv)])
                    soln_cv = clip.execute(pyclipper.INTERSECTION)
                    
                    #  Get the area of overlap with the DD chamber
                    A_CV += sum([pyclipper.area(loop) for loop in soln_cv])/scale_factor**2
                
                return A_CV
            
            t.append(theta)
            A.append(Atotal)
            Add.append(calculate_area(scaled_xdd, scaled_ydd))
            Ad1.append(calculate_area(scaled_xd1, scaled_yd1))
            Ac1_N.append(calculate_area(scaled_xc1_N, scaled_yc1_N))
            Ac1_Nm1.append(calculate_area(scaled_xc1_Nm1, scaled_yc1_Nm1))
            
            if plot:
                #ax.set_xlim(-0.025*scale_factor,0.025*scale_factor)
                #ax.set_ylim(-0.025*scale_factor,0.025*scale_factor)
                ax.set_aspect(1.0)  
                fig.savefig('disc_' + '{i:04d}'.format(i=i) + '.png')
        #  Save these values
        self.tdisc = np.array(t)
        self.Adisc_dd = np.array(Add)
        self.Adisc_d1 = np.array(Ad1)
        self.Adisc_c1_N = np.array(Ac1_N)
        self.Adisc_c1_Nm1 = np.array(Ac1_Nm1)
        
        if plot:
            fig = plt.figure()  
            ax = fig.add_subplot(111)
            ax.plot(t, (self.Adisc_dd+self.Adisc_d1+self.Adisc_c1_N+self.Adisc_c1_Nm1)*1e6, label='A')
            ax.plot(t, self.Adisc_dd*1e6, label='Add')
            ax.plot(t, self.Adisc_d1*1e6, label='Ad1')
            ax.plot(t, self.Adisc_c1_N*1e6, label='Ac1.N')
            ax.plot(t, self.Adisc_c1_Nm1*1e6, label='Ac1.(N-1)')
            ax.axvline(self.theta_d)
            #plt.legend()
            plt.xlabel('Crank angle [rad]')
            plt.ylabel('Area [mm$^2$]')
            fig.savefig('A_v_t.png')
            plt.close()
            
        if plot:
            print('making animation in disc_ani.gif')
            subprocess.check_call('convert disc_0*.png disc_ani.gif',shell=True)
            print('removing disc_0*.png')
            for file in glob.glob('disc_0*.png'):
                os.remove(file)
        
        #  Create spline interpolator objects
        self.spline_Adisc_DD = interp.splrep(t, self.Adisc_dd, k = 2, s = 0) # DD and port
        self.spline_Adisc_D1 = interp.splrep(t, self.Adisc_d1, k = 2, s = 0) #  D1 and port
        self.spline_Adisc_C1_N = interp.splrep(t, self.Adisc_c1_N, k = 2, s = 0) # C1_N and port
        self.spline_Adisc_C1_Nm1 = interp.splrep(t, self.Adisc_c1_Nm1, k = 2, s = 0) #C1_{N-1} and port
        
        print('done')
            
    @property
    def theta_d(self):
        return scroll_geo.theta_d(self.geo)
    
    @property
    def Vdisp(self):
        return -2*pi*self.geo.h*self.geo.rb*self.geo.ro*(3*pi
                                                         -2*self.geo.phi_fie
                                                         +self.geo.phi_fi0
                                                         +self.geo.phi_oo0)
    
    @property
    def Vratio(self):
        return ((3*pi-2*self.geo.phi_fie+self.geo.phi_fi0+self.geo.phi_oo0)
                /(-2*self.geo.phi_oos-3*pi+self.geo.phi_fi0+self.geo.phi_oo0))
    
    def V_injection(self, theta, V_tube = None):
        """
        Volume code for injection tube
        
        The tube volume can either be given by the keyword argument V_tube 
        (so you can easily have more than one injection tube), or it can be 
        provided by setting the Scroll class attribute V_inj_tube 
        and NOT providing the V_tube argument
        
        The injection tube volume is assumed to be constant, hence the derivative of volume is zero 
        """
        if V_tube is None:
            return self.V_inj_tube, 0.0
        else:
            return V_tube, 0.0
        
    def V_sa(self, theta, full_output=False):
        """
        Wrapper around the Cython code for sa calcs
        
        Parameters
        ----------
        theta: float
             angle in range [0,2*pi]
        
        Returns
        -------
        
        """
        return scroll_geo.SA(theta,self.geo,Vremove = self.geo.Vremove)[0:2]
        
    def V_s1(self,theta):
        """
        Wrapper around the Cython code for Vs1_calcs
        
        theta: angle in range [0,2*pi]
        """
        return scroll_geo.S1(theta,self.geo)[0:2]
        
    def V_s2(self,theta):
        """
        Wrapper around the Cython code for Vs1_calcs
        
        theta: angle in range [0,2*pi]
        """

        return scroll_geo.S2(theta,self.geo)[0:2]
    
    def V_c1(self,theta,alpha=1,full_output=False):
        """
        Wrapper around the Cython code for C1
        
        theta: angle in range [0,2*pi]
        alpha: index of compression chamber pair; 1 is for outermost set
        """
        return scroll_geo.C1(theta,alpha,self.geo)[0:2]
        
    def V_c2(self,theta,alpha=1,full_output=False):
        """
        Wrapper around the Cython code for C2
        
        theta: angle in range [0,2*pi]
        alpha: index of compression chamber pair; 1 is for outermost set
        """
        return scroll_geo.C2(theta,alpha,self.geo)[0:2]
        
    def V_d1(self,theta,full_output=False):
        """
        Wrapper around the compiled code for D1
        
        theta: angle in range [0,2*pi]
        """
        
        return scroll_geo.D1(theta,self.geo)[0:2]
    
    def V_d2(self,theta,full_output=False):
        """
        Wrapper around the compiled code for D2
        
        theta: angle in range [0,2*pi]
        """

        return scroll_geo.D2(theta,self.geo)[0:2]
    
    def V_dd(self,theta,full_output=False):
        """
        Wrapper around the compiled code for DD
        
        theta: angle in range [0,2*pi]
        alpha: index of compression chamber pair; 1 is for outermost set
        """
        if full_output==True:
            HTangles = {'1_i':None,'2_i':None,'1_o':None,'2_o':None}
            return scroll_geo.DD(theta,self.geo)[0:2],HTangles
        else:
            return scroll_geo.DD(theta,self.geo)[0:2]
        
    def V_ddd(self,theta,alpha=1,full_output=False):
        """
        Wrapper around the compiled code for DDD
        
        theta: angle in range [0,2*pi]
        alpha: index of compression chamber pair; 1 is for outermost set
        """
        if full_output==True:
            HTangles = {'1_i':None,'2_i':None,'1_o':None,'2_o':None}
            return scroll_geo.DDD(theta,self.geo)[0:2],HTangles
        else:
            return scroll_geo.DDD(theta,self.geo)[0:2]        
        
    def set_scroll_geo(self, *args, **kwargs):
        """
        Delegates to scroll_geo.set_scroll_geo to calculate the scroll geometry
        """

        set_scroll_geo(*args, geo=self.geo, **kwargs)
        
        #Set the flags to ensure all parameters are fresh
        self.__Setscroll_geo__=True
        self.__SetDiscGeo__=False
        
    def add_sensor(self, x, y):
        """
        Add a virtual sensor at the coordinates x,y
        
        Parameters
        ----------
        x : float
            Cartesian coordinates corresponding to the point on the fixed scroll
        y : float
            Cartesian coordinates corresponding to the point on the fixed scroll
        """
        if not hasattr(self.sensors,'coords'):
            self.sensors.coords = []
            self.sensors.T = []
            self.sensors.p = []
            self.sensors.rho = []
            
        self.sensors.coords.append((x,y))
        
    def determine_partner_CVs(self,x,y,N = 1000, theta = None):
        """
        For a given point, determine which control volume is connected to it 
        over the course of a rotation.  This can be useful for "instrumenting"
        of the numerical model
        
        Parameters
        ----------
        x : float
            X coordinate of the point
        y : float
            Y coordinate of the point
        N : int
            Number of elements in the crank angle array
        theta : iterable, optional
            The crank angles to be considered (N ignored if theta provided)
            
        Returns
        ------
        theta : numpy aray
            Crank angle array
        partners : list
            List of string keys for the CV found, ``None`` if no partner
        
        """
        if theta is not None:
            thetav = theta
        else:
            thetav = np.linspace(0,2*pi,1000)        
        
        #  The clipper library operates on integers, so we need to take our 
        #  floating point values and convert it to a large integer
        scale_factor = 1000000000
        
        # Scale to integers (a box one unit a side)
        scaled_x = [x*scale_factor,x*scale_factor+1,x*scale_factor+1,x*scale_factor,x*scale_factor]
        scaled_y = [y*scale_factor,y*scale_factor,y*scale_factor+1,y*scale_factor+1,y*scale_factor]
        
        from PDSim.misc.clipper import pyclipper
        
        partners = []
        for i, theta in enumerate(thetav):
            _found = False
            #  Find all the CVs that do have some overlap with this point
            for CVkey in self.CVs.keys:
                
                try:
                    # Get the coordinates of this chamber
                    xcv, ycv = scroll_geo.CVcoords(CVkey, self.geo, theta)
                    
                    # Get the coordinates in scaled values
                    scaled_xcv = xcv*scale_factor
                    scaled_ycv = ycv*scale_factor
                    
                    # Clip them
                    clip = pyclipper.Pyclipper()
                    clip.subject_polygon([pair for pair in zip(scaled_x, scaled_y)])
                    clip.clip_polygon([pair for pair in zip(scaled_xcv, scaled_ycv)])
                    
                    #  Actually do the intersection
                    soln = clip.execute(pyclipper.INTERSECTION)
                    
                    if soln:
                        partners.append(CVkey)
                        _found = True
                        break
                except:
#                    raise
                    pass
                
            if not _found:
                partners.append(None)
        
        return thetav, partners
    
    def poly_intersection_with_cvs(self, x, y, N, multiple_solns = 'sum', CVcoords = scroll_geo.CVcoords):
        """
        For a given polygon, determine the intersection area between a polygon
        and each of the control volumes in the compressor.  This can be useful
        to calculate the port open area for the discharge port, injection ports,
        pressure measurement taps, etc.
        
        Parameters
        ----------
        x : numpy-array-like
            X coordinates of the points
        y : numpy-array-like
            Y coordinates of the points
        N : int
            Number of elements in the crank angle array
        multiple_solns : str, optional, one of 'sum','error','warning', default is sum
            What to do when there are multiple intersection polygons
        CVcoords : function 
            The function that can provide the x,y coordinates for a given control volume,
            by default `scroll_geo.CVcoords`
            
        Returns
        ------
        theta : numpy aray
            Crank angle array
        area_dict : dictionary
            Dictionary with keys of keys of control volumes that have some 
            intersection, values are areas at each crank angle in ``theta`` 
            
        """
        
        #  The dictionary to store the overlap data
        area_dict = {}
        
        #  The crank angle array
        thetav = np.linspace(0, 2*np.pi, N)
        
        #  The clipper library operates on integers, so we need to take our 
        #  floating point values and convert it to a large integer
        scale_factor = 1000000000
        
        # Scale to integers
        scaled_xpoly = x*scale_factor
        scaled_ypoly = y*scale_factor
        
        #  Find all the CVs that do have some overlap with this polygon
        for CVkey in self.CVs.keys:
            
            #  Just try once to see if the key is acceptable
            try:
                xcv, ycv = CVcoords(CVkey, self.geo, 0)
            except KeyError:
                #  Not acceptable, skip this control volume
                continue
            
            Av = np.zeros_like(thetav)
            
            for i, theta in enumerate(thetav):
                
                #  Calculate the free area between the polygon and the chamber
                try:
                    xcv, ycv = scroll_geo.CVcoords(CVkey, self.geo, theta)
                except ValueError:
                    Av[i] = 0
                    continue
                    
                scaled_xcv = xcv*scale_factor
                scaled_ycv = ycv*scale_factor
                
                from PDSim.misc.clipper import pyclipper
                clip = pyclipper.Pyclipper()
                clip.subject_polygon([pair for pair in zip(scaled_xpoly, scaled_ypoly)])
                clip.clip_polygon([pair for pair in zip(scaled_xcv, scaled_ycv)])
                
                #  Actually do the intersection
                soln = clip.execute(pyclipper.INTERSECTION)
                
                #  Check the number of regions returned
                if len(soln) > 1:
                    if multiple_solns == 'error':
                        raise ValueError('More than one intersection region obtained in poly_intersection_with_cvs')
                    else:
                        msg = 'More than one intersection region obtained in poly_intersection_with_cvs'
                        warnings.warn(msg, warnings.RuntimeWarning)
                
                #  Sum up the solution regions
                A = 0
                for loop in soln:
                    scaled_x, scaled_y = zip(*loop)
                    x = [_/scale_factor for _ in scaled_x]
                    y = [_/scale_factor for _ in scaled_y]
                    A += scroll_geo.polyarea(x, y)
                    
                Av[i] = A
            
            #  If at least one value is non-zero, save an entry in the dictionary
            if np.sum(Av) > 0:
                area_dict[CVkey] = Av
        
        return thetav, area_dict
                    
    def set_disc_geo(self,Type,r2=0.0,**kwargs):
        """
        Set the discharge geometry for the scrolls
        
        Parameters
        ----------
        Type
            The type of 
        """
        if self.__Setscroll_geo__==False:
            raise ValueError("You must determine scroll wrap geometry by calling Setscroll_geo before setting discharge geometry.")
        
        #Use the compiled code
        scroll_geo.setDiscGeo(self.geo,Type,r2,**kwargs)
        
    def auto_add_CVs(self,inletState,outletState):
        """
        Adds all the control volumes for the scroll compressor.
        
        Parameters
        ----------
        inletState
            A :class:`State <CoolProp.State.State>` instance for the inlet to the scroll set.  Can be approximate
        outletState
            A :class:`State <CoolProp.State.State>` instance for the outlet to the scroll set.  Can be approximate
            
        Notes
        -----
        Uses the indices of 
        
        ============= ===================================================================
        CV            Description
        ============= ===================================================================
        ``sa``        Suction Area
        ``s1``        Suction chamber on side 1
        ``s2``        Suction chamber on side 2
        ``d1``        Discharge chamber on side 1
        ``d2``        Discharge chamber on side 2
        ``dd``        Central discharge chamber
        ``ddd``       Merged discharge chamber
        ``c1.i``      The i-th compression chamber on side 1 (i=1 for outermost chamber)
        ``c2.i``      The i-th compression chamber on side 2 (i=1 for outermost chamber)
        ============= ===================================================================
        """
        
        # Maximum number of pairs of compression chambers in existence
        nCmax = scroll_geo.nC_Max(self.geo)
        
        #Add all the control volumes that are easy.  Suction area and suction chambera
        self.add_CV(ControlVolume(key='sa',initialState=inletState.copy(),
                                VdVFcn=self.V_sa,becomes=['sa','s1','s2']))
        if nCmax > 0:
            self.add_CV(ControlVolume(key='s1',initialState=inletState.copy(),
                                VdVFcn=self.V_s1,becomes='c1.1'))
            self.add_CV(ControlVolume(key='s2',initialState=inletState.copy(),
                                VdVFcn=self.V_s2,becomes='c2.1'))
        else:
            self.add_CV(ControlVolume(key='s1',initialState=inletState.copy(),
                                VdVFcn=self.V_s1,becomes='s1'))
            self.add_CV(ControlVolume(key='s2',initialState=inletState.copy(),
                                VdVFcn=self.V_s2,becomes='s2'))
        
        # Discharge chambers are also easy.  Assume that you start with 'ddd' chamber merged.
        # No problem if this isn't true.
        self.add_CV(ControlVolume(key='d1',initialState=outletState.copy(),
                                VdVFcn=self.V_d1,exists=False))
        self.add_CV(ControlVolume(key='d2',initialState=outletState.copy(),
                                VdVFcn=self.V_d2,exists=False))
        self.add_CV(ControlVolume(key='dd',initialState=outletState.copy(),
                                VdVFcn=self.V_dd,exists=False))
        self.add_CV(ControlVolume(key='ddd',initialState=outletState.copy(),
                                VdVFcn=self.V_ddd,discharge_becomes='dd'))

        #Add each pair of compression chambers
        
        # Must have at least one pair
        #~ if nCmax < 1:
            #~ raise AssertionError('nCmax ['+str(nCmax)+'] must be at least 1')
        for alpha in range(1,nCmax+1):
            keyc1 = 'c1.'+str(alpha)
            keyc2 = 'c2.'+str(alpha)
            if alpha==1:
                #It is the outermost pair of compression chambers
                initState = State.State(inletState.Fluid,
                                        dict(T=inletState.T,
                                             D=inletState.rho)
                                        )
                
            else:
                #It is not the first CV, more involved analysis
                #Assume isentropic compression from the inlet state at the end of the suction process
                T1 = inletState.T
                s1 = inletState.s
                rho1 = inletState.rho
                k = inletState.cp/inletState.cv
                V1 = self.V_s1(2*pi)[0]
                V2 = self.V_c1(0,alpha)[0]
                #Mass is constant, so rho1*V1 = rho2*V2
                rho2 = rho1 * V1 / V2
                
                T2guess = T1*(V1/V2)**(k-1)
                # Now don't know temperature or pressure, but you can assume
                # it is isentropic to find the temperature
                temp = inletState.copy()
                def resid(T):
                    temp.update(dict(T=T, D=rho2))
                    return temp.s-s1
                # temp has now been updated
                initState=temp.copy()
            if alpha<nCmax:
                # Does not change definition at discharge angle
                disc_becomes_c1 = 'c1.'+str(alpha)
                disc_becomes_c2 = 'c2.'+str(alpha)
                # It is not the innermost pair of chambers, becomes another 
                # set of compression chambers at the end of the rotation
                becomes_c1 = 'c1.'+str(alpha+1)
                becomes_c2 = 'c2.'+str(alpha+1)
            else:
                #It is the innermost pair of chambers, becomes discharge chamber
                #at the discharge angle
                disc_becomes_c1 = 'd1'
                disc_becomes_c2 = 'd2'
                becomes_c1 = 'c1.'+str(alpha+1) #Not used - CV dies at disc.
                becomes_c2 = 'c2.'+str(alpha+1) #Not used - CV dies at disc.
                
            self.add_CV(ControlVolume(key=keyc1,
                                      initialState=initState.copy(),
                                      VdVFcn=self.V_c1,
                                      VdVFcn_kwargs={'alpha':alpha},
                                      discharge_becomes=disc_becomes_c1,
                                      becomes=becomes_c1))
            
            self.add_CV(ControlVolume(key=keyc2,
                                      initialState=initState.copy(),
                                      VdVFcn=self.V_c2,
                                      VdVFcn_kwargs={'alpha':alpha},
                                      discharge_becomes=disc_becomes_c2,
                                      becomes=becomes_c2))
    
    def auto_add_leakage(self, flankFunc = None, radialFunc = None, radialFunc_kwargs = {}, flankFunc_kwargs = {}):
        """
        Add all the leakage terms for the compressor
        
        Parameters
        ----------
        flankFunc : function
            The function to be used for the flank leakage path
        flankFunc_kwargs : function
            Dictionary of terms to be passed to the flank leakage function
        radialFunc : function
            The function to be used for the radial leakage path
        radialFunc_kwargs : dict
            Dictionary of terms to be passed to the radial leakage function
            
        """
        
        if flankFunc is not None:
            #Do the flank leakages
            self.auto_add_flank_leakage(flankFunc, flankFunc_kwargs)
            
        if radialFunc is not None:
            #Do the radial leakages
            self.auto_add_radial_leakage(radialFunc, radialFunc_kwargs)
        
    def auto_add_radial_leakage(self, radialFunc, radialFunc_kwargs):
        """
        A function to add all the radial leakage terms
        
        Parameters
        ----------
        radialFunc : function
            The function that will be called for each radial leakage
        radialFunc_kwargs : dict
            Dictionary of terms to be passed to the radial leakage function
        """
        #Get all the radial leakage pairs
        pairs = scroll_geo.radial_leakage_pairs(self.geo)
        
        # Loop over all the radial leakage pairs possible for the given geometry
        for pair in pairs:
            if ('sa' in pair or 's1' in pair or 's2' in pair) and hasattr(self,'disable_radial_suction') and self.disable_radial_suction:
                warnings.warn('radial s1-c1 disabled',RuntimeWarning)
                continue
            self.add_flow(FlowPath(key1=pair[0],
                                   key2=pair[1],
                                   MdotFcn=radialFunc,
                                   MdotFcn_kwargs = radialFunc_kwargs
                                   )
                          )
        
    def auto_add_flank_leakage(self, flankFunc, flankFunc_kwargs = {}):
        """
        A function to add all the flank leakage terms
        
        Parameters
        ----------
        flankFunc : function
            The function that will be called for each flank leakage
        flankFunc_kwargs : function
            Dictionary of terms to be passed to the flank leakage function
        """
        
        # Always a s1-c1 leakage and s2-c2 leakage
        if hasattr(self,'disable_flank_suction') and self.disable_flank_suction:
            warnings.warn('flank s1-c1 disabled',RuntimeWarning)
        else:
            self.add_flow(FlowPath(key1='s1',key2='c1.1',MdotFcn=flankFunc, MdotFcn_kwargs = flankFunc_kwargs))
            self.add_flow(FlowPath(key1='s2',key2='c2.1',MdotFcn=flankFunc, MdotFcn_kwargs = flankFunc_kwargs))
        
        # Only add the DDD-S1 and DDD-S2 flow path if there is one set of
        # compression chambers.   
        if scroll_geo.nC_Max(self.geo) == 1:
            self.add_flow(FlowPath(key1 = 's1', key2 = 'ddd',MdotFcn = self.DDD_to_S, MdotFcn_kwargs = flankFunc_kwargs))
            self.add_flow(FlowPath(key1 = 's2', key2 = 'ddd',MdotFcn = self.DDD_to_S, MdotFcn_kwargs = flankFunc_kwargs))
        
        #Add each pair of compression chambers
        nCmax = scroll_geo.nC_Max(self.geo)
        
        # Must have at least one pair
        # assert (nCmax>=1)
        
        for alpha in range(1, nCmax+1):
            keyc1 = 'c1.'+str(alpha)
            keyc2 = 'c2.'+str(alpha)
            
            if alpha <= nCmax - 1:
                #Leakage between compression chambers along a path
                self.add_flow(FlowPath(key1 = keyc1,
                                       key2 = 'c1.'+str(alpha+1),
                                       MdotFcn = flankFunc,
                                       MdotFcn_kwargs = flankFunc_kwargs))
                self.add_flow(FlowPath(key1 = keyc2,
                                       key2 = 'c2.'+str(alpha+1),
                                       MdotFcn = flankFunc, 
                                       MdotFcn_kwargs = flankFunc_kwargs))
                
            elif alpha==nCmax:
                #Leakage between the discharge region and the innermost chamber
                self.add_flow(FlowPath(key1 = keyc1, key2='ddd', MdotFcn=flankFunc, MdotFcn_kwargs = flankFunc_kwargs))
                self.add_flow(FlowPath(key1 = keyc2, key2='ddd', MdotFcn=flankFunc, MdotFcn_kwargs = flankFunc_kwargs))
            
            flankFunc_kwargs_copy = copy.deepcopy(flankFunc_kwargs)
            # Update the flag so that this term will only be evaluated when the number of pairs of 
            # compression chambers in existence will be equal to     
            flankFunc_kwargs_copy['Ncv_check'] = nCmax - 1
            
            if alpha == nCmax - 1:
                # Leakage between the discharge region and the next-most inner chamber when the innermost chambers
                # have been swallowed into the discharge region
                self.add_flow(FlowPath(key1 = keyc1, key2 = 'ddd', MdotFcn = flankFunc, MdotFcn_kwargs = flankFunc_kwargs_copy))
                self.add_flow(FlowPath(key1 = keyc2, key2 = 'ddd', MdotFcn = flankFunc, MdotFcn_kwargs = flankFunc_kwargs_copy))
                self.add_flow(FlowPath(key1 = keyc1, key2 = 'd1', MdotFcn = flankFunc, MdotFcn_kwargs = flankFunc_kwargs_copy))
                self.add_flow(FlowPath(key1 = keyc2, key2 = 'd2', MdotFcn = flankFunc, MdotFcn_kwargs = flankFunc_kwargs_copy))
    
    def calculate_scroll_mass(self):
        """
        Calculate the mass and centroid of the orbiting scroll
        
        Returns
        -------
        mtotal : float
            Total mass of the orbiting scroll (including any additional mass added at centroid)
            
        zcm__thrust_surface : float
            The location of the centroid above the thrust surface. z = 0 is at the thrust surface, and positive z direction is towards the orbiting scroll 
        """
        tplate = getattr(self.mech,'scroll_plate_thickness')
        rho = getattr(self.mech,'scroll_density')
        mplus = getattr(self.mech,'scroll_added_mass')
        Lbearing = getattr(self.mech,'L_crank_bearing',0)
        Dijournal = getattr(self.mech,'D_crank_bearing',0)
        Dplate = getattr(self.mech,'scroll_plate_diameter')
        Dojournal = 1.5*Dijournal
        
        Vwrap,cx,cy = scroll_geo.scroll_wrap(self.geo)
        
        mwrap = rho * Vwrap
        mplate = rho * pi * tplate * Dplate**2/4.0
        mjournal = rho * pi * Lbearing * (Dojournal**2-Dijournal**2)/4.0
        mtotal = mwrap + mplate + mjournal + mplus
        
        zwrap = tplate + self.geo.h/2
        zplate = tplate/2
        zjournal = -Lbearing/2
        
        zcm__thrust_surface = (mwrap*zwrap + mjournal*zjournal + mplate*zplate)/mtotal
        
        return mtotal, zcm__thrust_surface

    def heat_transfer_coefficient(self, key):
        
#        Pr=Pr_mix(Ref,Liq,T_avg,p_avg,xL_avg); //[-]
#        Re=4.0*mdot/2.0/(PI*mu_mix(Ref,Liq,T_avg,p_avg,xL_avg)*Dh); //[-]
#        hc=0.023*k_mix(Ref,Liq,T_avg,p_avg,xL_avg)/Dh*pow(Re,0.8)*pow(Pr,0.4); //[kW/m^2-K]
#        // Jang and Jeong correction for spiral geometry
#        f=scroll->States.omega/(2*PI);
#        Amax=scroll->geo.ro;
#        Ubar=scroll->massFlow.mdot_tot/(4*scroll->geo.ro*scroll->geo.hs*rho);
#        St=f*Amax/Ubar;
#        hc*=1.0+8.48*(1-exp(-5.35*St));
#        // Tagri and Jayaraman correction for transverse oscillation
#        r_c=scroll->geo.rb*(0.5*phi_1_i+0.5*phi_2_i-scroll->geo.phi.phi_fi0);
#        hc*=1.0+1.77*Dh/r_c;
        return 1.0

    def wrap_heat_transfer(self, **kwargs):
        """
        This function evaluates the anti-derivative of the differential of wall heat transfer, and returns the amount of scroll-wall heat transfer in kW
        
        Parameters
        ----------
        hc : float
            Heat transfer coefficient [kW/m2/K]
        hs : float
            Scroll wrap height [m]
        rb : float
            Base circle radius [m]
        phi1 : float
            Larger involute angle [rad]
        phi2 : float
            Smaller involute angle [rad]
        phi0 : float
            Initial involute angle [rad]
        T_scroll : float
            Lump temperature of the scroll wrap [K]
        T_CV : float
            Temperature of the gas in the CV [K]
        dT_dphi : float
            Derivative of the temperature along the scroll wrap [K/rad]
        phim : float
            Mean involute angle of wrap used for heat transfer [rad]
        
        Notes
        -----
        ``phi1`` and ``phi2`` are defined such that ``phi1`` is always the
        larger involute angle in value
        """
        #Use the compiled version from the cython code
        return _Scroll.involute_heat_transfer(self,**kwargs)
    
    def heat_transfer_callback(self, theta):
        """
        The scroll simulation heat transfer callback for HT to the fluid in the 
        chambers
        
        ``heat_transfer_callback`` for ``PDSimCore.derivs`` must be of the 
        form::
        
            heat_transfer_callback(theta)
            
        but we need to get the inlet and outlet states to get the linear 
        temperature profile in the scroll wrap. Thus we wrap the callback 
        we would like to call in this function that allows us to determine
        the inlet and outlet state at run-time.
        """
        State_inlet = self.Tubes.Nodes[self.key_inlet]
        State_outlet = self.Tubes.Nodes[self.key_outlet]
        return self._heat_transfer_callback(theta, State_inlet, State_outlet)
    
    def _heat_transfer_callback(self, theta, State_inlet, State_outlet, HTC_tune = 1.0):
        """
        A private function to actually do the heat transfer analysis
        """
        # dT_dphi is generally negative because as you move to the 
        # outside of the scroll (larger phi), the temperature goes down because
        # you are moving towards low pressure and low temperature

        Tsuction = State_inlet.T
        Tdischarge = State_outlet.T
        dT_dphi = (Tsuction - Tdischarge) / (self.geo.phi_fie - self.geo.phi_oos)
        phim = 0.5*self.geo.phi_fie + 0.5*self.geo.phi_oos
        
        Q = []
        for key in self.CVs.exists_keys:
            Q.append(self.calcHT(theta,key,HTC_tune,dT_dphi,phim))
        return arraym(Q)
        
    def HT_angles(self, theta, geo, key):
        return symm_scroll_geo.HT_angles(theta, self.geo, key)
        
    def calcHT(self, theta, key, HTC_tune, dT_dphi, phim): 
        
        #TODO: calculate HTC
        hc = self.HTC #[kW/m2/K]
        
        ## If HT is turned off, no heat transfer
        if abs(hc) < 1e-10 or HTC_tune <= 0.0 or key.startswith('inj') or key == 'sa' or key == 'dd':
            return 0.0
        elif key == 'ddd':
            # ddd is a combination of the heat transfer in the d1, d2, and
            # dd chambers
            Q_d1 = self.calcHT(theta,str('d1'),HTC_tune,dT_dphi,phim)
            Q_d2 = self.calcHT(theta,str('d2'),HTC_tune,dT_dphi,phim)
            return Q_d1 + Q_d2
                
        
            
        #Get the bounding angles for the control volume
        angles = self.HT_angles(theta, self.geo, key)
        
        if angles is None:
            return 0.0
        
        T_scroll = self.Tlumps[0]
        T_CV = self.CVs[key].State.T
        # The heat transfer rate of the inner involute on 
        # the outer wrap of the chamber
        Q_outer_wrap = scroll_geo.involute_heat_transfer(hc, 
                                                   self.geo.h, 
                                                   self.geo.rb, 
                                                   angles.phi_1_i, 
                                                   angles.phi_2_i, 
                                                   angles.phi_i0, 
                                                   T_scroll,
                                                   T_CV, 
                                                   dT_dphi, 
                                                   phim)
        
        # The heat transfer rate of the outer involute on 
        # the inner wrap of the chamber
        Q_inner_wrap = scroll_geo.involute_heat_transfer(hc, 
                                                   self.geo.h, 
                                                   self.geo.rb, 
                                                   angles.phi_1_o, 
                                                   angles.phi_2_o, 
                                                   angles.phi_o0,
                                                   T_scroll,
                                                   T_CV, 
                                                   dT_dphi, 
                                                   phim)
        
        return HTC_tune *(Q_outer_wrap + Q_inner_wrap)
        
    def step_callback(self,t,h,Itheta):
        """
        Here we test whether the control volumes need to be
        a) Merged
        b) Adjusted because you are at the discharge angle
        
        """ 
        # This gets called at every step, or partial step
        self.theta = t
        
        def angle_difference(angle1,angle2):
            # Due to the periodicity of angles, you need to handle the case where the
            # angles wrap around - suppose theta_d is 6.28 and you are at an angles of 0.1 rad
            #, the difference should be around 0.1, not -6.27
            # 
            # This brilliant method is from http://blog.lexique-du-net.com/index.php?post/Calculate-the-real-difference-between-two-angles-keeping-the-sign
            # and the comment of user tk
            return (angle1-angle2+pi)%(2*pi)-pi
        
        def IsAtMerge(eps = 0.001, eps_d1_higher = 0.01, eps_dd_higher = 0.00001):
            pressures = [self.CVs['d1'].State.p,
                         self.CVs['d2'].State.p,
                         self.CVs['dd'].State.p]
            p_max = max(pressures)
            p_min = min(pressures)
            if abs(p_min/p_max-1)<eps_dd_higher:
                return True
            # For over compression cases, the derivatives don't tend to drive
            # the pressures together, and as a result you need to relax the 
            # convergence quite a bit
            elif angle_difference(t, scroll_geo.theta_d(self.geo))>1.2 and abs(p_min/p_max-1)<eps_d1_higher:
                return True
            else:
                return False
            
        disable=False
        
        if t<self.theta_d<t+h and self.__before_discharge__==False:
            # Take a step almost up to the discharge angle
            #print 'almost to discharge'
            disable = True
            h = self.theta_d-t-1e-15
            self.__before_discharge__ = True
            
        elif self.__before_discharge__ == True:
            # At the discharge angle
            #print 'At the discharge angle'
            ########################
            #Reassign chambers
            ########################
            #Find chambers with a discharge_becomes flag
            for key in self.CVs.exists_keys:
                if self.CVs[key].discharge_becomes in self.CVs.keys:
                    #Set the state of the "new" chamber to be the old chamber
                    oldCV=self.CVs[key]
                    if not oldCV.exists==True: raise AttributeError("old CV doesn't exist")
                    
                    newCV=self.CVs[oldCV.discharge_becomes]
                    newCV.State.update({'T':oldCV.State.T,'D':oldCV.State.rho})
                    oldCV.exists=False
                    newCV.exists=True
            
            self.__before_discharge__=False
            
            self.update_existence()
            
            #Re-calculate the CV volumes
            V,dV = self.CVs.volumes(self.theta_d+1e-10)
            #Update the matrices using the new CV definitions
            self.T[self.CVs.exists_indices,Itheta] = self.CVs.T
            self.p[self.CVs.exists_indices,Itheta] = self.CVs.p
            self.m[self.CVs.exists_indices,Itheta] = arraym(self.CVs.rho)*V
            self.rho[self.CVs.exists_indices,Itheta] = arraym(self.CVs.rho)
            
            # This has to be quite large - why?
            h = 2e-10
            disable='no_integrate' # This means that no actual update will be made, simply a copy from the old CV to the new CV
       
        elif self.CVs['d1'].exists and IsAtMerge():
            
            #Build the volume vector using the old set of control volumes (pre-merge)
            V,dV=self.CVs.volumes(t)
            
            print('merging')
            
            if self.__hasLiquid__==False:

                #Density
                rhod1=self.CVs['d1'].State.rho
                rhod2=self.CVs['d2'].State.rho
                rhodd=self.CVs['dd'].State.rho
                #Density
                pd1=self.CVs['d1'].State.p
                pd2=self.CVs['d2'].State.p
                pdd=self.CVs['dd'].State.p
                #Internal energy
                ud1=self.CVs['d1'].State.u
                ud2=self.CVs['d2'].State.u
                udd=self.CVs['dd'].State.u
                #Internal energy
                Td1=self.CVs['d1'].State.T
                Td2=self.CVs['d2'].State.T
                Tdd=self.CVs['dd'].State.T
                #Volumes
                Vdict=dict(zip(self.CVs.exists_keys,V))
                Vd1=Vdict['d1']
                Vd2=Vdict['d2']
                Vdd=Vdict['dd']
                
                Vddd=Vd1+Vd2+Vdd
                m=rhod1*Vd1+rhod2*Vd2+rhodd*Vdd
                U_before=ud1*rhod1*Vd1+ud2*rhod2*Vd2+udd*rhodd*Vdd
                rhoddd=m/Vddd
                #guess the mixed temperature as a volume-weighted average
                T=(Td1*Vd1+Td2*Vd2+Tdd*Vdd)/Vddd
                p=(pd1*Vd1+pd2*Vd2+pdd*Vdd)/Vddd
                #Must conserve mass and internal energy (instantaneous mixing process)
                Fluid = self.CVs['ddd'].State.Fluid
                
                temp = self.CVs['ddd'].State.copy()
                def resid(T):
                    temp.update(dict(T=T,D=rhoddd))
                    return temp.u - U_before/m
                    
                T_u = optimize.newton(resid, T)
                
                self.CVs['ddd'].State.update({'T':T_u,'D':rhoddd})
                U_after=self.CVs['ddd'].State.u*self.CVs['ddd'].State.rho*Vddd
                
                DeltaU=m*(U_before-U_after)
                if abs(DeltaU)>1e-5:
                    raise ValueError('Internal energy not sufficiently conserved in merging process')
                
                self.CVs['d1'].exists=False
                self.CVs['d2'].exists=False
                self.CVs['dd'].exists=False
                self.CVs['ddd'].exists=True
                
                self.update_existence()
                
                #Re-calculate the CV
                V,dV=self.CVs.volumes(t)
                self.T[self.CVs.exists_indices,Itheta] = self.CVs.T
                self.p[self.CVs.exists_indices,Itheta] = self.CVs.p
                self.m[self.CVs.exists_indices,Itheta] = arraym(self.CVs.rho)*V
                self.rho[self.CVs.exists_indices,Itheta] = arraym(self.CVs.rho)
                
            else:
                raise NotImplementedError('no flooding yet')
            disable=True 
              
        elif t>self.theta_d:
            self.__before_discharge__=False
            disable=False
            
        return disable,h
        
    def crank_bearing(self, W):
        
        JB = journal_bearing(r_b = self.mech.D_crank_bearing/2,
                             L = self.mech.L_crank_bearing,
                             omega = self.omega,
                             W = W,
                             c = self.mech.c_crank_bearing,
                             eta_0 = self.mech.mu_oil
                             )
        self.losses.crank_bearing_dict = JB
    
        return JB['Wdot_loss']/1000.0
        
    def upper_bearing(self, W):
        """
        Moment balance around the upper bearing gives the force for
        the lower bearing.  Torques need to balance around the upper bearing
        """
        
        JB = journal_bearing(r_b = self.mech.D_upper_bearing/2,
                             L = self.mech.L_upper_bearing,
                             omega = self.omega,
                             W = W,
                             c = self.mech.c_upper_bearing,
                             eta_0 = self.mech.mu_oil
                             )
        self.losses.upper_bearing_dict = JB

        return JB['Wdot_loss']/1000.0
    
    def lower_bearing(self, W):
        """
        Moment balance around the upper bearing gives the force for
        the lower bearing.  Torques need to balance around the upper bearing
        """
        
        JB = journal_bearing(r_b = self.mech.D_lower_bearing/2,
                             L = self.mech.L_lower_bearing,
                             omega = self.omega,
                             W = W,
                             c = self.mech.c_lower_bearing,
                             eta_0 = self.mech.mu_oil
                             )
        self.losses.lower_bearing_dict = JB

        return JB['Wdot_loss']/1000.0
    
    def thrust_bearing(self):
        """
        The thrust bearing analysis
        """
        from PDSim.core.bearings import thrust_bearing
        V = self.geo.ro*self.omega
        #Use the corrected force to account for the decrease in back area due to the bearing
        N = self.forces.summed_Fz*1000 #[N]
        TB = thrust_bearing(mu = self.mech.thrust_friction_coefficient,
                            V = V,
                            N = N)
        self.losses.thrust_bearing_dict = TB
        return TB['Wdot_loss']/1000.0
    
    def mechanical_losses(self, shell_pressure = 'low:shell'):
        """
        Calculate the mechanical losses in the bearings
        
        Parameters
        ----------
            shell_pressure : string, 'low', 'low:shell', 'mid', or 'high'
            
            low uses the upstream pressure of the machine,
            
            low:shell uses the pressure after the inlet tube
            
            mid uses the average of upstream and downstream pressures
            
            high uses the pressure downstream of the machine 
            
        Returns
        -------
            Wdot_losses: float
            
                The total mechanical losses in kW

        """
        
        #inlet pressure [kPa]
        inlet_pressure = self.Tubes.Nodes[self.key_inlet].p
        outlet_pressure = self.Tubes.Nodes[self.key_outlet].p
        
        # Find the tube with the inlet node
        Tube = self.Tubes[self.key_inlet]
        # Get the state that is not the inlet state
        if Tube.key1 == b'self.key_inlet':
            shell_pressure_val = Tube.State2.p
        else:
            shell_pressure_val = Tube.State1.p
        
        # Get the shell pressure based on either the inlet or outlet pressure
        # based on whether it is a low-pressure or high-pressure shell
        if shell_pressure == 'low':
            back_pressure = min((inlet_pressure, outlet_pressure))
        elif shell_pressure == 'low:shell':
            back_pressure = min((shell_pressure_val, outlet_pressure))
        elif shell_pressure == 'high':
            back_pressure = max((inlet_pressure, outlet_pressure))
        elif shell_pressure == 'mid':
            back_pressure = (inlet_pressure + outlet_pressure)/2
        else:
            raise KeyError("keyword argument shell_pressure must be one of 'low', 'low:shell', 'mid' or 'high'; received '"+str(shell_pressure)+"'")
        
        #Calculate the force terms: force profiles, mean values, etc. 
        self.calculate_force_terms(orbiting_back_pressure = back_pressure)
        
        if not hasattr(self,'losses'):
            self.losses = struct()
            
        if not hasattr(self.mech,'journal_tune_factor'):
            warnings.warn('mech.journal_tune_factor not provided; using 1.0')
            self.mech.journal_tune_factor = 1.0
            
        if hasattr(self.mech,'specified_mechanical_efficiency') and isinstance(self.mech.specified_mechanical_efficiency, float):
            return (1-self.mech.specified_mechanical_efficiency)*self.Wdot_pv
        elif hasattr(self.mech,'specified_mechanical_losses_kW') and isinstance(self.mech.specified_mechanical_losses_kW, float):
            return self.mech.specified_mechanical_losses_kW
        elif hasattr(self.mech,'detailed_analysis') and self.mech.detailed_analysis == True:
            self.detailed_mechanical_analysis()
            return self.forces.Wdot_total_mean # [kW]
        else:
            #Conduct the calculations for the bearings [N]
            W_OSB = np.sqrt((self.forces.summed_Fr + self.forces.inertial)**2+self.forces.summed_Ft**2)*1000
            
            _slice = list(range(self.Itheta+1))
            
            self.losses.crank_bearing = np.zeros_like(W_OSB)
            self.losses.upper_bearing = np.zeros_like(W_OSB)
            self.losses.lower_bearing = np.zeros_like(W_OSB)
            
            for i in _slice:
                self.losses.crank_bearing[i] = self.crank_bearing(W_OSB[i])
                self.losses.upper_bearing[i] = self.upper_bearing(W_OSB[i]*(1+1/self.mech.L_ratio_bearings))
                self.losses.lower_bearing[i] = self.lower_bearing(W_OSB[i]/self.mech.L_ratio_bearings)
                
            self.losses.thrust_bearing = self.thrust_bearing()
            
            theta = self.t[_slice]
            theta_range = theta[-1]-theta[0]
            
            # Sum up each loss v. theta curve
            self.losses.summed = self.losses.crank_bearing + self.losses.upper_bearing + self.losses.lower_bearing + self.losses.thrust_bearing
            
            # Get the mean losses over one cycle
            self.losses.bearings  = np.trapz(self.losses.summed[_slice], theta)/theta_range
            
            print('mechanical losses: ', self.losses.bearings)
            return self.losses.bearings #[kW]
    
    def post_cycle(self):
        # Run the base-class method to set HT terms, etc. - also calls lumps_energy_balance_callback
        PDSimCore.post_cycle(self)
        
        # Update the heat transfer to the gas in the shell
        self.suction_heating()
    
    def post_solve(self):
        
        # Run the base-class method to set HT terms, etc. - also calls lumps_energy_balance_callback
        PDSimCore.post_solve(self)
        
        # Build the pressure profiles
        self.build_pressure_profile()
        
        # Build the virtual sensor profile
        self.build_sensor_profile()
        
        # Add some more entries to the summary
        self.summary.eta_oi = self.eta_oi
        self.summary.Wdot_electrical = self.Wdot_electrical
        
    def build_sensor_profile(self):
        """
        Build the thermo data for each point along the process
        """
        
        # First check if this run uses any virtual sensors
        if not hasattr(self.sensors,'T') and not hasattr(self.sensors,'coords'): return
        
        print('preparing to calculate sensor profiles, this could take a while')
        for x,y in self.sensors.coords:
            theta, partners = self.determine_partner_CVs(x, y, theta = self.t)

            ##  Collect all the data            
            T,p,rho = [],[],[]
            for i,partner in enumerate(partners):
                if partner is None:
                    T.append(np.nan)
                    p.append(np.nan)
                    rho.append(np.nan)
                else:
                    ##  Get integer key for partner
                    I = self.CVs.index(partner)
                    
                    T.append(self.T[I,i])
                    p.append(self.p[I,i])
                    rho.append(self.rho[I,i])
            
            self.sensors.T.append(np.array(T))
            self.sensors.p.append(np.array(p))
            self.sensors.rho.append(np.array(rho))
        
    def Nc_max(self):
        return [scroll_geo.nC_Max(self.geo) for i in range(2)]
        
    
    def build_pressure_profile(self):
        """
        Build the pressure profile, tracking along s1,c1.x,d1,ddd and
        s2,c2.x,d2,ddd and store them in the variables summary.theta_profile,
        summary.p1_profile, summary.p2_profile
        """
        
        # Calculate along one path to track one set of pockets through the whole process
        theta = self.t
        pcopy = self.p.copy() # make a copy to avoid risk of pointer-to-data-overwriting problem

        # Suction chambers
        p1 = pcopy[self.CVs.index('s1')]
        p2 = pcopy[self.CVs.index('s2')]
        
        Nc_max1, Nc_max2 = self.Nc_max()
        
        assert len(theta) == len(p1) == len(p2)
        
        for path, Nc_max in zip([1,2],[Nc_max1, Nc_max2]):
            if Nc_max > 1:
                for alpha in range(1,Nc_max):
                    # Compression chambers up to the next-to-innermost set are handled
                    # just like the suction chambers
                    if path == 1:
                        theta = np.append(theta, self.t + 2*pi*alpha)
                        p1 = np.append(p1, pcopy[self.CVs.index('c1.'+str(alpha))])
                    else:
                        p2 = np.append(p2, pcopy[self.CVs.index('c2.'+str(alpha))])
                    
            # Innermost compression chamber begins to be tricky
            # By definition innermost compression chamber doesn't make it to the 
            # end of the rotation
            next_theta = self.t + 2*pi*Nc_max
            if path == 1:
                next_p1 = pcopy[self.CVs.index('c1.'+str(Nc_max))]
                next_p1[np.isnan(next_p1)] = 0
            else:
                next_p2 = pcopy[self.CVs.index('c2.'+str(Nc_max))]
                next_p2[np.isnan(next_p2)] = 0
        
        pd1 = pcopy[self.CVs.index('d1')]
        pd2 = pcopy[self.CVs.index('d2')]
        pddd = pcopy[self.CVs.index('ddd')]
        
        assert len(theta) == len(p1) == len(p2)
        
        # Now check if d1 and d2 end before the end of the rotation (they don't 
        # neccessarily)
        if np.isnan(pd1[0]) and np.isnan(pd1[self.Itheta]):
            # d1 & d2 end before the end of the rotation
            # straightforward analysis (just add on pd1 and pd2)
            pd1[np.isnan(pd1)] = 0
            pd2[np.isnan(pd2)] = 0
            next_p1 += pd1
            next_p2 += pd2
            
            # So we know that ddd DOES exist at the beginning/end of the rotation
            # work backwards to find the first place that the ddd does exist
            pdddA = pddd.copy()
            pdddB = pddd.copy()
            
            i = self.Itheta
            while i > 0:
                if np.isnan(pdddA[i]):
                    i += 1
                    break;
                i -= 1
            pdddA[0:i] = 0 # This is the end of the rotation
            next_p1 += pdddA
            next_p2 += pdddA
            
            theta = np.append(theta, next_theta)
            p1 = np.append(p1, next_p1)
            p2 = np.append(p2, next_p2)
            
            i = 0
            while i < len(pdddB):
                if np.isnan(pdddB[i]):
                    break;
                i += 1
            
            pdddB[i::] = np.nan # This is the beginning of the next rotation
            
            theta = np.append(theta, self.t + 2*pi*(Nc_max1 + 1))
            p1 = np.append(p1, pdddB)
            p2 = np.append(p2, pdddB)
        
        # Now check if d1 & d2 still exist at the end of the rotation
        elif not np.isnan(pd1[0]) and not np.isnan(pd1[self.Itheta]):
            # d1 & d2 don't end before the end of the rotation
            pd1A = pd1.copy()
            pd1B = pd1.copy()
            pd2A = pd2.copy()
            pd2B = pd2.copy()
            
            i = self.Itheta
            while i > 0:
                if np.isnan(pd2A[i]):
                    i += 1
                    break;
                i -= 1
            pd1A[0:i] = 0 # This is the end of the rotation
            pd2A[0:i] = 0 # This is the end of the rotation
            next_p1 += pd1A
            next_p2 += pd2A
            
            theta = np.append(theta, next_theta)
            p1 = np.append(p1, next_p1)
            p2 = np.append(p2, next_p2)
            
            last_theta = self.t + 2*pi*(Nc_max1 + 1)
            last_p1 = pddd.copy()
            last_p2 = pddd.copy()
            last_p1[np.isnan(last_p1)] = 0
            last_p2[np.isnan(last_p2)] = 0
            
            i = 0
            while i < len(pd1B):
                if np.isnan(pd1B[i]):
                    break;
                i += 1
            if i == len(pd1B)-1:
                raise ValueError('d1B could not find NaN')
            
            pd1B[i::] = 0
            pd2B[i::] = 0
            last_p1 += pd1B
            last_p2 += pd2B
            
            theta = np.append(theta, last_theta)
            p1 = np.append(p1, last_p1)
            p2 = np.append(p2, last_p2)
        
        self.summary.theta_profile = theta
        self.summary.p1_profile = p1
        self.summary.p2_profile = p2
        
        assert len(theta) == len(p1) == len(p2)
        
    def ambient_heat_transfer(self, Tshell):
        """
        The amount of heat transfer from the compressor to the ambient [kW]
        """
        return self.h_shell*self.A_shell*(Tshell-self.Tamb)
    
    def initial_motor_losses(self, eta_a = 0.8):
        """
        Assume an adiabatic efficiency to estimate the motor power and motor losses
        """
        
        for Tube in self.Tubes:
            if self.key_inlet in [Tube.key1, Tube.key2]:
                mdot = Tube.mdot
                
        inletState = self.Tubes.Nodes[self.key_inlet]
        outletState = self.Tubes.Nodes[self.key_outlet]
        s1 = inletState.s
        h1 = inletState.h
        temp = inletState.copy()
        temp.update(dict(S = s1, P = outletState.p))
        h2s = temp.h
        
        if outletState.p > inletState.p:
            #Compressor Mode
            h2 = h1 + (h2s-h1)/eta_a
        else:
            #Expander Mode
            h2 = h1 + (h2s-h1)*eta_a
        
        # A guess for the compressor mechanical power based on 70% efficiency [kW]
        Wdot = abs(mdot*(h2-h1))
        
        if self.motor.type == 'const_eta_motor':
            eta = self.motor.eta_motor
        else:
            #The efficiency and speed [-,rad/s] from the mechanical power output
            eta, self.omega = self.motor.invert_map(Wdot)
        
        #Motor losses [kW]
        self.motor.losses = Wdot*(1.0/eta-1)
        
    def suction_heating(self):
        if hasattr(self,'motor'):
            # If some fraction of heat from motor losses is going to get added
            # to suction flow
            if 0.0 <= self.motor.suction_fraction <= 1.0:
                for Tube in self.Tubes:
                    # Find the tube that has one of the keys starting with 'inlet'
                    if Tube.key1.startswith('inlet') or Tube.key2.startswith('inlet'):
                        #Add some fraction of the motor losses to the inlet gas 
                        Tube.Q_add = self.motor.losses * self.motor.suction_fraction
                    else:
                        Tube.Q_add = 0.0
                        
    def pre_run(self, N = 40000):
        """
        Intercepts the call to pre_run and does some scroll processing, then 
        calls the base class function
        """
        
        #  Get an initial guess before running at all for the motor losses.
        self.initial_motor_losses() #set the parameter self.motor.losses
        
        #  Run the suction heating code
        self.suction_heating()
        
        #  Calculate the dischare port free area
        self.cache_discharge_port_blockage(xport = self.geo.xvec_disc_port, 
                                           yport = self.geo.yvec_disc_port)
        
        #  Set the index for each control volume 
        for FP in self.Flows:
            FP.key1Index = scroll_geo.get_compressor_CV_index(FP.key1)
            FP.key2Index = scroll_geo.get_compressor_CV_index(FP.key2)
        
        #  Call the base class function        
        PDSimCore.pre_run(self, N = N)
        
        self.CVs['sa'].ForceFcn = symm_scroll_geo.SA_forces
        self.CVs['s1'].ForceFcn = symm_scroll_geo.S1_forces
        self.CVs['s2'].ForceFcn = symm_scroll_geo.S2_forces
        self.CVs['d1'].ForceFcn = symm_scroll_geo.D1_forces
        self.CVs['d2'].ForceFcn = symm_scroll_geo.D2_forces
        self.CVs['dd'].ForceFcn = symm_scroll_geo.DD_forces
        self.CVs['ddd'].ForceFcn = symm_scroll_geo.DDD_forces
        for i in range(10):
            try:
                Fcn = lambda theta, geo: symm_scroll_geo.C1_forces(theta, i+1, geo)
                self.CVs['c1.' + str(i+1)].ForceFcn = Fcn
            except:
                pass
            try:
                Fcn = lambda theta, geo : symm_scroll_geo.C2_forces(theta, i+1, geo)
                self.CVs['c2.' + str(i+1)].ForceFcn = Fcn
            except:
                pass
        
    def guess_lump_temps(self, T0):
        """
        Guess the temperature of the lump
        
        Parameters
        ----------
        T0 : float
            First guess for temperature [K]
        """
        
        # First try to just alter the lump temperature with the gas heat transfer
        # rate fixed
        
        def OBJECTIVE(x):
            self.Tlumps[0] = x
            #Run the tubes
            for tube in self.Tubes:
                tube.TubeFcn(tube)
            #
            return self.lump_energy_balance_callback()[0]
        
        print(OBJECTIVE(T0-50))
        print(OBJECTIVE(T0+50))
        return optimize.newton(OBJECTIVE, T0)
        
    def lump_energy_balance_callback(self):
        """
        The callback where the energy balance is carried out on the lumps
        
        Notes
        -----
        Derivation for electrical power of motor:
        
        .. math ::
            
            \\eta _{motor} = \\frac{\\dot W_{shaft}}{\\dot W_{shaft} + \\dot W_{motor}}
            
        .. math ::
            
            {\\eta _{motor}}\\left( \\dot W_{shaft} + \\dot W_{motor} \\right) = \\dot W_{shaft}
            
        .. math::
        
            \\dot W_{motor} = \\frac{\\dot W_{shaft}}{\\eta _{motor}} - \\dot W_{shaft}
        """
        
        #For the single lump
        # HT terms are positive if heat transfer is TO the lump
        Qnet = 0.0
        Qnet -= sum([Tube.Q for Tube in self.Tubes])
        
        self.Qamb = self.ambient_heat_transfer(self.Tlumps[0])
        
        self.mech.Wdot_losses = self.mechanical_losses('low:shell') 
        # Shaft power from forces on the orbiting scroll from the gas in the pockets [kW]
        self.Wdot_forces = self.omega*self.forces.mean_tau
#        if self.geo.is_symmetric():
#            
#        else:
#            import warnings
#            warnings.warn('No ML for asymmetric for now')
#            self.mech.Wdot_losses = 0.0
        
        # Heat transfer with the ambient; Qamb is positive if heat is being removed, thus flip the sign
        Qnet -= self.Qamb
        
        Qnet += self.mech.Wdot_losses
        # Heat transfer with the gas in the working chambers.  mean_Q is positive
        # if heat is transfered to the gas in the working chamber, so flip the 
        # sign for the lump
        Qnet -= self.HTProcessed.mean_Q
        
        
        
        self.Wdot_mechanical = self.Wdot_pv + self.mech.Wdot_losses
        
        #The actual torque required to do the compression [N-m]
        self.tau_mechanical = self.Wdot_mechanical / self.omega * 1000
        
        # 2 Options for the motor losses:
        # a) Constant efficiency
        # b) Based on the torque-speed-efficiency motor
        
        if self.motor.type == 'const_eta_motor':
            self.eta_motor = self.motor.eta_motor
        elif self.motor.type == 'motor_map':
            # Use the motor map to calculate the slip rotational speed [rad/s]
            # and the motor efficiency as a function of the torque [N-m]
            eta, omega = self.motor.apply_map(self.tau_mechanical)
            self.eta_motor = eta
            self.omega = omega
        else:
            raise AttributeError('motor.type must be one of "const_eta_motor" or "motor_map"')
        
        #Motor losses [kW]
        self.motor.losses = self.Wdot_mechanical*(1/self.eta_motor-1)
        
        #Electrical Power
        self.Wdot_electrical = self.Wdot_mechanical + self.motor.losses
        
        if hasattr(self,'Wdot_i'):
            #Overall isentropic efficiency
            self.eta_oi = self.Wdot_i/self.Wdot_electrical
        
#        #Set the heat input to the suction line
#        self.suction_heating()
        
        if self.verbosity > 0:
            print('At this iteration')
            print('    Electrical power:', self.Wdot_electrical,'kW')
            print('    Mass flow rate:', self.mdot,'kg/s')
            if hasattr(self,'Wdot_i'):
                print('    Over. isentropic:', self.eta_oi,'-')
            if hasattr(self,'eta_v'):
                print('    Volumetric:', self.eta_v,'-')
        
        #Want to return a list
        return [Qnet]

    def TubeCode(self,Tube,**kwargs):
        Tube.Q = flow_models.IsothermalWallTube(Tube.mdot,
                                                Tube.State1,
                                                Tube.State2,
                                                Tube.fixed,
                                                Tube.L,
                                                Tube.ID,
                                                T_wall=self.Tlumps[0],
                                                Q_add = Tube.Q_add,
                                                alpha = Tube.alpha
                                                )
        
    
    def DDD_to_S(self,FlowPath,flankFunc = None,**kwargs):
        if  flankFunc is None:
            flankFunc = self.FlankLeakage
        # If there are any compression chambers, don't evaluate this flow
        # since the compression chambers "get in the way" of flow directly from 
        # ddd to s1 and s2
        if scroll_geo.getNc(self.theta,self.geo) > 0:
            return 0.0
        else:
            return flankFunc(FlowPath)
            
    def D_to_DD(self,FlowPath,X_d =1.0,**kwargs):
        
        FlowPath.A=scroll_geo.Area_d_dd(self.theta,self.geo)
        try:
            return flow_models.IsentropicNozzle(FlowPath.A,
                                                FlowPath.State_up,
                                                FlowPath.State_down)*X_d
        except ZeroDivisionError:
            return 0.0
        
    def DISC_DD(self, FP, X_d = 1.0):
        """
        The flow path function for the flow between discharge port and the DD chamber
        """
        
        FP.A = interp.splev(self.theta, self.spline_Adisc_DD)
        try:
            return flow_models.IsentropicNozzle(FP.A, FP.State_up, FP.State_down) * X_d
        except ZeroDivisionError:
            return 0.0
        
    def DISC_D1(self, FP, X_d = 1.0):
        """
        The flow path function for the flow between discharge port and the D1 chamber
        """
        
        FP.A = interp.splev(self.theta, self.spline_Adisc_D1)
        try:
            return flow_models.IsentropicNozzle(FP.A, FP.State_up, FP.State_down) * X_d
        except ZeroDivisionError:
            return 0.0
    
    def DISC_C1_N(self, FP, X_d = 1.0):
        """
        The flow path function for the flow between discharge port and the D1 chamber
        """
        
        FP.A = interp.splev(self.theta, self.spline_Adisc_C1_N)
        try:
            return flow_models.IsentropicNozzle(FP.A, FP.State_up, FP.State_down) * X_d
        except ZeroDivisionError:
            return 0.0
            
    def DISC_C1_Nm1(self, FP, X_d = 1.0):
        """
        The flow path function for the flow between discharge port and the D1 chamber
        """
        
        FP.A = interp.splev(self.theta, self.spline_Adisc_C1_Nm1)
        try:
            return flow_models.IsentropicNozzle(FP.A, FP.State_up, FP.State_down) * X_d
        except ZeroDivisionError:
            return 0.0
     
    def SA_S1(self, FlowPath, X_d = 1.0, X_d_precompression = 1.0):
        """
        A wrapper for the flow between the suction area and the S1 chamber
        
        Notes
        -----
        If geo.phi_ie_offset is greater than 0, the offset geometry will be 
        used to calculate the flow area.  Otherwise the conventional analysis 
        will be used.
        """
        V,dV = scroll_geo.S1(self.theta,self.geo)[0:2]
        
        if dV >= 0:
            #  Normal incoming suction flow
            _X_d = X_d 
        else:
            #  Back flow at the end
            _X_d = X_d_precompression
            
        if abs(self.geo.phi_ie_offset) > 1e-12:
            FlowPath.A = scroll_geo.Area_s_s1_offset(self.theta, self.geo)
        else:
            FlowPath.A = scroll_geo.Area_s_sa(self.theta, self.geo)
             
        try:
            mdot = _X_d*flow_models.IsentropicNozzle(FlowPath.A,
                                                FlowPath.State_up,
                                                FlowPath.State_down)
            return mdot
        except ZeroDivisionError:
            return 0.0   
        
    def SA_S2(self, *args, **kwargs):
        """
        A thin wrapper to the default suction area-suction flow
        """
        return self.SA_S(*args,**kwargs)
        
    def SA_S(self, FlowPath, X_d = 1.0, X_d_precompression = 1.0):
        
        V,dV = scroll_geo.S1(self.theta,self.geo)[0:2]
        
        if dV >= 0:
            #  Normal incoming suction flow
            _X_d = X_d 
        else:
            #  Back flow at the end
            _X_d = X_d_precompression
            
        FlowPath.A = _X_d*scroll_geo.Area_s_sa(self.theta, self.geo)
        try:
            mdot = flow_models.IsentropicNozzle(FlowPath.A,
                                                FlowPath.State_up,
                                                FlowPath.State_down)
            return mdot
        except ZeroDivisionError:
            return 0.0
        
#    def _get_injection_CVkey(self,phi,theta,inner_outer):
#        """
#        Find the CV that is in contact with the given injection port location
#        
#        Parameters
#        ----------
#        phi : float
#            Involute angle of the injection port location
#        theta : float
#            Crank angle in radians in the range [:math:`0,2\pi`]
#        inner_outer : string ['i','o']
#            'i' : involute angle corresponds to outer surface of fixed scroll
#            'o' : involute angle corresponds to inner surface of orb. scroll 
#            
#        Notes
#        -----
#        Typically 'i' will require a positive offset in involute angle of 
#        :math:`\pi` radians
#        """
#        if inner_outer == 'i':
#            phi_0 = self.geo.phi_i0
#            phi_s = self.geo.phi_is
#            phi_e = self.geo.phi_ie
#        elif inner_outer == 'o':
#            phi_0 = self.geo.phi_o0
#            phi_s = self.geo.phi_os
#            phi_e = self.geo.phi_oe-pi # The actual part of the wrap that could 
#                                       # have an injection port 
#        
#        Nc = scroll_geo.getNc(theta, self.geo)    
#        #Start at the outside of the given scroll wrap
#        # x1 where x is s,d,c has the inner involute of the fixed scroll as 
#        # its outer surface
#        if phi_e > phi > phi_e-theta:     
#            #It is a suction chamber    
#            return 's1' if inner_outer == 'i' else 's2'
#            
#        elif phi_e-theta > phi > phi_e-theta-2*pi*Nc:
#            #It is one of the compression chambers, figure out which one
#            for I in range(Nc+1):
#                if phi_e - theta - 2*pi*(I-1) > phi > phi_e - theta - 2*pi*I:
#                    i_str = '.'+str(I)
#                    break
#            return 'c1'+i_str if inner_outer == 'i' else 'c2'+i_str
#        
#        else:
#            return 'd1' if inner_outer == 'i' else 'd2'
#        
#    def Injection_to_Comp(self,FlowPath,phi,inner_outer,check_valve = False, A = 7e-6, **kwargs):
#        """
#        Function to calculate flow rate between injection line and chamber
#        
#        Parameters
#        ----------
#        FlowPath : FlowPath instance
#        phi : involute angle where the port is located
#        inner_outer : string ['i','o']
#            'i' : involute angle corresponds to outer surface of fixed scroll
#            'o' : involute angle corresponds to inner surface of orb. scroll 
#        check_valve : boolean
#            If ``True``, there is an idealized check valve and flow can only go 
#            from chambers with key names that start with `injCV` to other chambers.
#            If ``False``, flow can go either direction
#        
#        """
#        #1. Figure out what CV is connected to the port
#        partner_key = self._get_injection_CVkey(phi, self.theta, inner_outer)
#        FlowPath.A = A
#        #2. Based on what CV is connected to the port, maybe quit
#        if partner_key in ['d1', 'd2'] and 'ddd' in [FlowPath.key_up, 
#                                                     FlowPath.key_down]:
#            # Other chamber based on geometry is d1 or d2 but they are not 
#            # defined due to the angle but ddd is, and as a result, use 
#            # ddd
#            #
#            # Don't do anything, just let it go to the next section even though
#            # 'd1' or 'd2 is not key_up or key_down 
#            pass
#        
#        elif partner_key not in [FlowPath.key_up, FlowPath.key_down]:
#            return 0.0
#        # If the pressure in the injection line is below the other chamber and 
#        # you are using a theoretical check valve with instantaneous closing, 
#        # then there is no back flow, and hence no flow at all
#        elif check_valve:
#            if FlowPath.key_down.startswith('inj'):
#                return 0.0
#            
##                #This will be negative
##                DELTAp = FlowPath.State_down.p - FlowPath.State_up.p
##            else:
##                #This will be positive
##                DELTAp = FlowPath.State_up.p - FlowPath.State_down.p 
##            
##            # Using an approximation to a Heaviside step function to close off the
##            # port gradually and improve numerical convergence due to continuous
##            # first derivative
##            if -10 < DELTAp < 10.0:
##                FlowPath.A *=  1/(1+np.exp(-10*(DELTAp-2)))
##            elif DELTAp < -10.0:
##                return 0.0
#        
#        mdot = flow_models.IsentropicNozzle(FlowPath.A,
#                                            FlowPath.State_up,
#                                            FlowPath.State_down)
#        return mdot
        
    def calculate_force_terms(self,
                              orbiting_back_pressure=None):
        """
        Calculate the force profiles, mean forces, moments, etc.
        
        Parameters
        ----------
        orbiting_back_pressure : float, or class instance
            If a class instance, must provide a function __call__ that takes as its first input the Scroll class
        
        """
        if not isinstance(orbiting_back_pressure,float):
            raise ValueError('calculate_force_terms must get a float back pressure for now')
            
        if not hasattr(self.mech,'scroll_plate_thickness'):
            warnings.warn('"mech.scroll_plate_thickness" not found, using 2*scroll wrap thickness')
            self.mech.scroll_plate_thickness = 2*self.geo.t
        
        if not hasattr(self.mech,'scroll_zcm__thrust_surface'):
            warnings.warn('"mech.scroll_zcm__thrust_surface" not found, using 0')
            self.mech.scroll_zcm__thrust_surface = 0
        
        self.forces = struct()
        
        #Get the slice of indices that are in use.  At the end of the simulation
        #execution this will be the full range of the indices, but when used
        # at intermediate iterations it will be a subset of the indices
        _slice = list(range(self.Itheta+1))
        
        t = self.t[_slice]
        
        # Remove control volumes that have no change in volume as they are not "real" control volumes
        isrealCV_val = (np.max(self.V[:,_slice], axis = 1) - np.min(self.V[:,_slice], axis = 1))
        
        ####################################################
        #############  Inertial force components ###########
        ####################################################
        
        #: The magnitude of the inertial forces on the orbiting scroll [kN]
        self.forces.inertial = self.mech.orbiting_scroll_mass * self.omega**2 * self.geo.ro / 1000
        
        ####################################################
        #############  Normal force components #############
        ####################################################
        # The force of the gas in each chamber pushes the orbiting scroll away
        # from the working chambers
        # It is only the active slice
        self.forces.Fz = self.p[:,_slice]*self.V[:,_slice]/self.geo.h
        
        # Remove all the NAN placeholders and replace them with zero values
        self.forces.Fz[np.isnan(self.forces.Fz)] = 0
        
        for i, val in enumerate(isrealCV_val):
            if not np.isnan(val) and val < 1e-14:
                self.forces.Fz[i,:] = 0
        
        # Sum the terms for the applied gas force from each of the control volumes
        self.forces.summed_Fz = np.sum(self.forces.Fz, axis = 0) #kN
#        
        # The back gas pressure on the orbiting scroll pushes the scroll back down
        # Subtract the back pressure from all the control volumes 
        self.forces.Fbackpressure = orbiting_back_pressure*self.V[:,_slice]/self.geo.h
        
        # Remove all the NAN placeholders and replace them with zero values
        self.forces.Fbackpressure[np.isnan(self.forces.Fbackpressure)] = 0
        
        for i, val in enumerate(isrealCV_val):
            if not np.isnan(val) and val < 1e-14:
                self.forces.Fbackpressure[i,:] = 0
        
        self.forces.summed_Fbackpressure = np.sum(self.forces.Fbackpressure, axis = 0) #kN
        
        self.forces.summed_Fz -= self.forces.summed_Fbackpressure
        
        # Add the net axial force generated by the gas at the top of the scroll wrap
        self.forces.Faxial_involute = self.scroll_involute_axial_force(t) 
        
        self.forces.summed_Faxial_involute = np.sum(self.forces.Faxial_involute, axis = 0)
        
        self.forces.summed_Fz += self.forces.summed_Faxial_involute
        
        # Calculate the mean axial force
        self.forces.mean_Fz = np.trapz(self.forces.summed_Fz, t)/(2*pi)
        
        ####################################################
        #############  "Radial" force components ###########
        ####################################################
        self.forces.Fx = np.zeros((self.CVs.N,len(self.t[_slice])))
        self.forces.Fy = np.zeros_like(self.forces.Fx)
        self.forces.fxp = np.zeros_like(self.forces.Fx)
        self.forces.fyp = np.zeros_like(self.forces.Fx)
        self.forces.cx = np.zeros_like(self.forces.Fx)
        self.forces.cy = np.zeros_like(self.forces.Fx)
        self.forces.Mz = np.zeros_like(self.forces.Fx)
        
        # A map of CVkey to function to be called to get force components
        # All functions in this map use the same call signature and are "boring"
        # Each function returns a dictionary of terms

        for CVkey in self.CVs.keys:
            try:
                geo_components = [self.CVs[CVkey].ForceFcn(theta, self.geo) for theta in self.t[_slice]]
            except BaseException as BE:
                print('no forces for', CVkey, 'with error:', BE)             
                geo_components = []
                
            if geo_components:
                I = self.CVs.index(CVkey)
                p = self.p[I,_slice]
                self.forces.fxp[I,:] = [comp['fx_p'] for comp in geo_components]
                self.forces.fyp[I,:] = [comp['fy_p'] for comp in geo_components]
                self.forces.Fx[I,:] = [comp['fx_p'] for comp in geo_components]*p
                self.forces.Fy[I,:] = [comp['fy_p'] for comp in geo_components]*p
                self.forces.cx[I,:] = [comp['cx'] for comp in geo_components]
                self.forces.cy[I,:] = [comp['cy'] for comp in geo_components]
                self.forces.Mz[I,:] = [comp['M_O_p'] for comp in geo_components]*p
        
        # The thrust load from JUST the working chambers
        self.forces.summed_gas_Fz = np.sum(self.forces.Fz, axis = 0)
        
        # Point of action of the thrust forces - weighted by the axial force
        self.forces.cx_thrust = np.sum(self.forces.cx*self.forces.Fz, axis = 0) / self.forces.summed_gas_Fz
        self.forces.cy_thrust = np.sum(self.forces.cy*self.forces.Fz, axis = 0) / self.forces.summed_gas_Fz
        
        self.forces.THETA = self.geo.phi_fie-pi/2-self.t[_slice]
        # Position of the pin as a function of crank angle
        self.forces.xpin = self.geo.ro*np.cos(self.forces.THETA)
        self.forces.ypin = self.geo.ro*np.sin(self.forces.THETA)
        
        # Remove all the NAN placeholders
        self.forces.Fx[np.isnan(self.forces.Fx)] = 0
        self.forces.Fy[np.isnan(self.forces.Fy)] = 0
        self.forces.Fz[np.isnan(self.forces.Fz)] = 0
        self.forces.Mz[np.isnan(self.forces.Mz)] = 0
        
        # Sum the terms at each crank angle
        self.forces.summed_Fx = np.sum(self.forces.Fx, axis = 0) #kN
        self.forces.summed_Fy = np.sum(self.forces.Fy, axis = 0) #kN
        
        # Moment around the x axis and y-axis from the thrust load caused by the working chambers
        # Sign convention on force is backwards here (in this case, forces pointing down are positive
        # r is vector from (xpin,ypin) to (cx,cy), F is -Fz*k
        #        |   i          j      k   |  
        #   M =  | cx-xpin   cy-ypin   0   | = - (cy-ypin) * Fz *i + (cx- xpin) * Fz * j
        #        |   0          0     -Fz  |
        #   
        self.forces.Mx = -(self.forces.cy - self.forces.ypin)*self.forces.Fz
        self.forces.My = +(self.forces.cx - self.forces.xpin)*self.forces.Fz
        
        # Moment around the x-axis and y-axis from the applied gas load on the orbiting scroll wrap relative to the thrust plane
        #        |   i    j     k  |  
        #   M =  |   0    0     z  | = - Fy*z *i + Fx * z * j
        #        |   Fx   Fy    0  |
        #
        self.forces.Mx += -self.forces.Fy*(self.mech.scroll_plate_thickness + self.geo.h/2)
        self.forces.My += self.forces.Fx*(self.mech.scroll_plate_thickness + self.geo.h/2)
        
        self.forces.summed_Mz = np.sum(self.forces.Mz, axis = 0) #kN-m
        self.forces.summed_Mx = np.sum(self.forces.Mx, axis = 0) #kN-m
        self.forces.summed_My = np.sum(self.forces.My, axis = 0) #kN-m
        
        # Magnitude of the overturning moment generated by the gas forces (Fr, Ftheta, Fz)
        self.forces.Moverturn_gas = np.sqrt(self.forces.summed_Mx**2+self.forces.summed_My**2)
        
        # The moments from the backpressure acting on the back side of the orbiting scroll
        # They must be added on separately because otherwise they are added in NCV times,
        # which is not the proper behavior
        # r is vector from (xpin,ypin) to (0,0), F is Fbp*k
        #        |   i        j      k   |  
        #   M =  | 0-xpin   0-ypin   0   | = -ypin * Fbp *i + xpin * Fbp * j
        #        |   0        0     Fbp  |
        #        
        self.forces.summed_Mx += -self.forces.ypin*self.forces.summed_Fbackpressure
        self.forces.summed_My +=  self.forces.xpin*self.forces.summed_Fbackpressure

        # Moment around the x-axis and y-axis from the inertial force of the orbiting scroll
        # They must be added on separately because otherwise they are added in NCV times,
        # which is not the proper behavior
        #
        # Components of the inertial force in the x- and y-directions [kN]
        Fcx = self.forces.inertial*np.cos(self.forces.THETA)
        Fcy = self.forces.inertial*np.sin(self.forces.THETA)
        
        # Inertial overturning moment acts through the center of mass of the orbiting scroll
        # Moment around the x-axis and y-axis from the applied gas load on the orbiting scroll 
        # wrap relative to the thrust plane
        #
        #        |   i      j     k  |  
        #   M =  |   0      0     z  | = - Fcy*z *i + Fcx * z * j
        #        |   Fcx   Fcy    0  |
        #
        self.forces.summed_Mx += -Fcy*(self.mech.scroll_zcm__thrust_surface)
        self.forces.summed_My += Fcx*(self.mech.scroll_zcm__thrust_surface)
        
        # Moment around the x-axis and y-axis from the journal bearing applied forces on the orbiting scroll
        self.forces.Fx_bearing_simple = -(np.cos(self.forces.THETA)*self.forces.inertial + self.forces.summed_Fx)
        self.forces.Fy_bearing_simple = -(np.sin(self.forces.THETA)*self.forces.inertial + self.forces.summed_Fy)
        #        |   i      j     k   |  
        #   M =  |   0      0     -h  | = Fy*h *i - Fx * h * j
        #        |   Fx     Fy    0   |
        #
        if hasattr(self.mech,'D_crank_bearing'):
            self.forces.summed_Mx +=  self.forces.Fy_bearing_simple*(self.mech.D_crank_bearing/2)
            self.forces.summed_My += -self.forces.Fx_bearing_simple*(self.mech.D_crank_bearing/2)
        else:
            self.forces.summed_Mx += 0
            self.forces.summed_My += 0
        
        self.forces.Moverturn_thrust = np.sum(self.forces.Fz*np.sqrt((self.forces.cy - self.forces.ypin)**2+(self.forces.cx - self.forces.xpin)**2),axis=0)
        self.forces.Moverturn_gas = (self.mech.scroll_plate_thickness + self.geo.h/2)*np.sum(np.sqrt(self.forces.Fx**2+self.forces.Fy**2),axis=0)
        self.forces.Moverturn_inertia = (self.mech.scroll_zcm__thrust_surface)*np.sqrt(Fcy**2+Fcx**2)
        if hasattr(self.mech,'D_crank_bearing'):
            self.forces.Moverturn_bearing = (self.mech.D_crank_bearing/2)*np.sqrt(self.forces.Fx_bearing_simple**2+self.forces.Fy_bearing_simple**2)
        else:
            self.forces.Moverturn_bearing = None
        self.forces.Moverturn = np.sqrt(self.forces.summed_Mx**2+self.forces.summed_My**2)
        
#        import matplotlib.pyplot as plt
#        plt.plot(t,self.forces.Moverturn_thrust,'b')
#        plt.plot(t,self.forces.Moverturn_gas,'r')
#        plt.plot(t,self.forces.Moverturn_inertia,'g')
#        plt.plot(t,self.forces.Moverturn_bearing,'k')
#        plt.show()
        
        # Center of reaction in the global coordinate system
        self.forces.x_thrust_reaction = -self.forces.summed_My/self.forces.summed_Fz+self.forces.xpin #Fz is positive if pushing UP on the OS
        self.forces.y_thrust_reaction = -self.forces.summed_Mx/self.forces.summed_Fz+self.forces.ypin #Fz is positive if pushing UP on the OS
        self.forces.r_thrust_reaction = np.sqrt(self.forces.x_thrust_reaction**2+ self.forces.y_thrust_reaction**2)
        
        #Calculate the radial force on the crank pin at each crank angle
        #The radial component magnitude is just the projection of the force onto a vector going from origin to center of orbiting scroll
        self.forces.Fr = (np.cos(self.forces.THETA)*self.forces.Fx + np.sin(self.forces.THETA)*self.forces.Fy)
        #
        #Components of the unit vector in the direction of rotation
        x_dot = +np.sin(self.forces.THETA)
        y_dot = -np.cos(self.forces.THETA)
        # Direction of rotation is opposite the positive theta direction, so need to flip sign for Ft
        self.forces.Ft = -(x_dot*self.forces.Fx+y_dot*self.forces.Fy)
        
        #Remove all the NAN placeholders
        self.forces.Fr[np.isnan(self.forces.Fr)] = 0
        self.forces.Ft[np.isnan(self.forces.Ft)] = 0
        #Sum the terms at each crank angle
        self.forces.summed_Fr = np.sum(self.forces.Fr,axis = 0) #kN
        self.forces.summed_Ft = np.sum(self.forces.Ft,axis = 0) #kN

        self.forces.Fm = np.sqrt(self.forces.summed_Fx**2+self.forces.summed_Fy**2)
        
        self.forces.tau = self.forces.xpin*self.forces.summed_Fy-self.forces.ypin*self.forces.summed_Fx
        # Calculate the mean normal force on the crank pin
        # This assumes a quasi-steady bearing where the film is well-behaved
        self.forces.mean_Fm = np.trapz(self.forces.Fm, self.t[_slice])/(2*pi)
        self.forces.mean_Fr = np.trapz(self.forces.summed_Fr, self.t[_slice])/(2*pi)
        self.forces.mean_Ft = np.trapz(self.forces.summed_Ft, self.t[_slice])/(2*pi)
        self.forces.mean_tau = np.trapz(self.forces.tau, self.t[_slice])/(2*pi)
        self.forces.mean_Mz = np.trapz(self.forces.summed_Mz, self.t[_slice])/(2*pi)
        
    def detailed_mechanical_analysis(self):
        """
        In this function the detailed mechanical analsysis is carried out.
        
        Notes
        -----
        This function implements the method of the documentation on the orbiting scroll
        forces
        
        The applied force on each of the bearings can be given from::
        
                        **|  |**    -----
                        **|  |**     |
                      |      |       | z_upper
                      |      |       |
                    **|      |**     |
                    **|      |**    -----   
                    **|      |**     |
                      |      |       |
                      |      |       |
                      |      |       |
                      |      |       | z_lower
                      |      |       |
                      |      |       |
                      |      |       |
                    **|      |**     |
                    **|      |**   -----
        
        and the ratio of bearing distances is given by ``mech.L_ratio_bearings`` which is z_lower/z_upper
        """
        
        if not hasattr(self,'losses'):
            self.losses = struct()
            
        muthrust = self.mech.thrust_friction_coefficient
        beta = self.mech.oldham_rotation_beta
        mu1 = mu2 = mu3 = mu4 = mu5 = self.mech.oldham_key_friction_coefficient
        r1 = r2 = r3 = r4 = self.mech.oldham_ring_radius
        w1 = w2 = w3 = w4 = self.mech.oldham_key_width
        mOR = self.mech.oldham_mass
        mOS = self.mech.orbiting_scroll_mass
        wOR = self.mech.oldham_thickness
        hkeyOR = self.mech.oldham_key_height
        
        # Fill in terms if they are not provided for backwards compatability
        for term in ['pin1_ybeta_offset','pin2_ybeta_offset','pin3_xbeta_offset','pin4_xbeta_offset']:
            if not hasattr(self.mech, term):
                warnings.warn('"mech.'+term+'" not found, using 0')
                setattr(self.mech,term,0)

        F1_ybeta_offset = self.mech.pin1_ybeta_offset
        F2_ybeta_offset = self.mech.pin2_ybeta_offset
        F3_xbeta_offset = self.mech.pin3_xbeta_offset
        F4_xbeta_offset = self.mech.pin4_xbeta_offset
        
        # Gravitional acceleration
        g = 9.80665 
        
        _slice = list(range(self.Itheta+1))
        theta = self.t[_slice]
        
        # The initial guesses for the moment generated by the journal bearing - 
        # it should be positive since Ms is negative and Ms and M_B act in 
        # opposite directions
        self.forces.F_B0 = np.sqrt((self.forces.summed_Fr + self.forces.inertial)**2+self.forces.summed_Ft**2)
        self.forces.mu_B = np.zeros_like(self.forces.F_B0)
        self.forces.mu_Bupper = np.zeros_like(self.forces.F_B0)
        self.forces.mu_Blower = np.zeros_like(self.forces.F_B0)
        self.forces.M_B0 = self.forces.mu_B*self.mech.D_crank_bearing/2*self.forces.F_B0
        
        THETA = self.geo.phi_fie-pi/2-theta
        vOR_xbeta = self.geo.ro*self.omega*(np.sin(THETA)*np.cos(beta)-np.cos(THETA)*np.sin(beta)) #Velocity of the oldham ring in the xbeta direction
        aOR_xbeta = self.geo.ro*self.omega**2*(-np.cos(THETA)*np.cos(beta)-np.sin(THETA)*np.sin(beta))
        UPSILON = vOR_xbeta/np.abs(vOR_xbeta)
        vOS_ybeta = -self.geo.ro*self.omega*(np.sin(THETA)*np.sin(beta)+np.cos(THETA)*np.cos(beta)) #Velocity of the orbiting scroll in the y_beta direction
        PSI = vOS_ybeta/np.abs(vOS_ybeta)
        vOS = self.geo.ro*self.omega
        aOS_x = -self.geo.ro*self.omega**2*np.cos(THETA)
        aOS_y = -self.geo.ro*self.omega**2*np.sin(THETA)
        
        Nsteps = self.Itheta+1
        A = np.zeros((4,4,Nsteps))
        b = np.zeros((4,Nsteps))
        self.forces.Fkey = np.zeros((4,Nsteps))
        
        # Make a matrix stack where each entry in the third index corresponds to a 4x4 matrix of the terms
        # Oldham x-beta direction
        A[0,0,:] = -mu1*UPSILON
        A[0,1,:] = -mu2*UPSILON
        A[0,2,:] = 1
        A[0,3,:] = -1
        b[0,:] = mOR*aOR_xbeta/1000
        
        # Oldham ybeta direction
        A[1,0,:] = 1    
        A[1,1,:] = -1 
        A[1,2,:] = -mu3*PSI
        A[1,3,:] = -mu4*PSI
        b[1,:] = 0
            
        # Oldham moments around the central z-direction axis
        A[2,0,:] = r1-mu1*UPSILON*(w1/2-F1_ybeta_offset)
        A[2,1,:] = r2+mu2*UPSILON*(w2/2+F2_ybeta_offset)
        A[2,2,:] = -r3+mu3*PSI*(w3/2-F3_xbeta_offset)
        A[2,3,:] = -r4-mu4*PSI*(w4/2+F4_xbeta_offset)
        b[2,:] = 0
        
        # Orbiting scroll moments around the central axis
        A[3,0,:] = 0
        A[3,0,:] = 0
        A[3,0,:] = r3-mu3*PSI*(w3/2-F3_xbeta_offset)
        A[3,0,:] = r4+mu4*PSI*(w4/2+F4_xbeta_offset)
        
        # Use the initial guess for the bearing moments and applied force
        self.forces.M_B = self.forces.M_B0
        self.forces.F_B = self.forces.F_B0
        
        step = 1
        
        # In the first step, we use the guess values for the bearing normal force, 
        # in the second step we use the back-calculated bearing normal force 
        while step <= 2:
            
            #Calculate the friction coefficient for each bearing
            for i in _slice:
                self.crank_bearing(self.forces.F_B[i]*1000)
                self.forces.mu_B[i] = self.losses.crank_bearing_dict['f']
                self.upper_bearing(self.forces.F_B[i]*1000*(1+1/self.mech.L_ratio_bearings))
                self.forces.mu_Bupper[i] = self.losses.upper_bearing_dict['f']
                self.lower_bearing(self.forces.F_B[i]*1000/self.mech.L_ratio_bearings)
                self.forces.mu_Blower[i] = self.losses.lower_bearing_dict['f']
                
            self.forces.M_B = self.forces.mu_B*self.mech.D_crank_bearing/2*self.forces.F_B
            
            # This term depends on M_B which is re-calculated at each iteration.  All other terms are independent of M_B
            b[3,:] = -self.forces.summed_Mz-self.forces.M_B-muthrust*(self.forces.summed_My*np.cos(THETA)-self.forces.summed_Mx*np.sin(THETA))
            
            # Walk through the stack of arrays, and obtain F1,F2,F3,F4 for each
            # crank angle
            for i in _slice:
                self.forces.Fkey[:,i] = np.linalg.solve(A[:,:,i], b[:,i])
        
                # from PDSim.plot.plots import debug_plots
                # debug_plots(self)
            F1, F2, F3, F4 = self.forces.Fkey
        
            # Bearing forces on the scroll re-calculated based on force balances in the x- and y-axes
            Fbx = mOS*aOS_x/1000-muthrust*self.forces.summed_Fz*np.sin(THETA)+mu3*PSI*F3*np.sin(beta)+mu4*PSI*F4*np.sin(beta)-F4*np.cos(beta)+F3*np.cos(beta)-self.forces.summed_Fx
            Fby = mOS*aOS_y/1000-muthrust*self.forces.summed_Fz*np.cos(THETA)-mu3*PSI*F3*np.cos(beta)-mu4*PSI*F4*np.cos(beta)-F4*np.sin(beta)+F3*np.sin(beta)-self.forces.summed_Fy
            Fbold = np.sqrt(Fbx**2+Fby**2)
            self.forces.M_B = self.forces.mu_B*self.mech.D_crank_bearing/2*Fbold
            self.forces.M_Bupper = self.forces.mu_Bupper*self.mech.D_upper_bearing/2*Fbold*(1.0+1.0/self.mech.L_ratio_bearings)
            self.forces.M_Blower = self.forces.mu_Blower*self.mech.D_lower_bearing/2*Fbold/self.mech.L_ratio_bearings
            
            self.forces.F_B = Fbold
            
            step += 1
        
        self.forces.Wdot_F1 = np.abs(F1*vOR_xbeta*mu1)
        self.forces.Wdot_F2 = np.abs(F2*vOR_xbeta*mu2)
        self.forces.Wdot_F3 = np.abs(F3*vOS_ybeta*mu3)
        self.forces.Wdot_F4 = np.abs(F4*vOS_ybeta*mu4)
        self.forces.Wdot_OS_journal = np.abs(self.omega*self.forces.M_B)*self.mech.journal_tune_factor
        self.forces.Wdot_upper_journal = np.abs(self.omega*self.forces.M_Bupper)*self.mech.journal_tune_factor
        self.forces.Wdot_lower_journal = np.abs(self.omega*self.forces.M_Blower)*self.mech.journal_tune_factor
        
        self.forces.Wdot_thrust = np.abs(muthrust*self.forces.summed_Fz*vOS)
        
        self.forces.Wdot_total = (  self.forces.Wdot_F1
                                  + self.forces.Wdot_F2
                                  + self.forces.Wdot_F3
                                  + self.forces.Wdot_F4
                                  + self.forces.Wdot_OS_journal
                                  + self.forces.Wdot_upper_journal 
                                  + self.forces.Wdot_lower_journal 
                                  + self.forces.Wdot_thrust)
        
        self.forces.Wdot_total_mean = np.trapz(self.forces.Wdot_total, theta)/(2*pi)
        print(self.forces.Wdot_total_mean,'average mechanical losses')
            
        import matplotlib.pyplot as plt
            
#        fig = plt.figure()
#        fig.add_subplot(111)
#        plt.plot(theta,self.forces.Wdot_F1,label='F1')
#        plt.plot(theta,self.forces.Wdot_F2,label='F2')
#        plt.plot(theta,self.forces.Wdot_F3,label='F3')
#        plt.plot(theta,self.forces.Wdot_F4,label='F4')
#        plt.plot(theta,self.forces.Wdot_OS_journal,label='journal')
#        plt.plot(theta,self.forces.Wdot_thrust,label='thrust')
#        plt.plot(theta,self.forces.Wdot_total,lw=3)
#        plt.legend()
#        plt.show()
#        
#        plt.close('all')
#        fig = plt.figure()
#        fig.add_subplot(141)
#        plt.gca().set_title('Oldham acceleration [g]')
#        plt.plot(theta,(F[2,:].T-F[3,:].T)*1000/9.81/mOR,'o',lw=2)
#        plt.plot(theta, aOR_xbeta/9.81,'-',lw=2)
#        fig.add_subplot(142)
#        plt.plot(theta,F.T)
#        plt.gca().set_title('Key forces [kN]')
#        fig.add_subplot(143)
#        plt.plot(theta,self.forces.summed_Ft)
#        plt.gca().set_title('Tangential force [kN]')
#        fig.add_subplot(144)
#        plt.plot(theta,self.forces.summed_Fr)
#        plt.gca().set_title('Radial force [kN]')
#        plt.show()
        
    def scroll_involute_axial_force(self, theta, p_backpressure = 0):
        """
        Calculate the axial force generated by the pressure distribution 
        along the top of the scroll wrap.  The force profile returned is the NET 
        force obtained by subtracting the back pressure from the applied force
        
        Pressure along inner and outer walls is considered to act over one-half 
        of the thickness of the scroll wrap.
        
        Notes
        -----
        The following assumptions are employed:
        
        1. Involute extended to the base circle to account for discharge region
        2. Half of the width of the scroll is assumed to see the upstream pressure and the other half sees the downstream pressure
        
        The length of an involute section can be given by
        
        .. math::
            
            s = r_b\\left(\\frac{\\phi_2^2-\\phi_1^2}{2}-\\phi_{0}(\\phi_2-\\phi_1)\\right)
        
        Returns
        -------
        F: numpy array
            Axial force matrix from the top of the scroll generated by each control volume [kN]
        """
        
        _slice = list(range(len(theta)))
        
        # Get the break angle (simplified solution)
        phi_s_sa = self.geo.phi_ooe-pi
        
        # Get the number of compression chambers in existence at each crank angle
        nC = np.array([scroll_geo.getNc(t,self.geo) for t in theta])
        
        F = np.zeros_like(self.p)
        F = F[:,_slice]
        
        # Parameters for the SA chamber
        phi2 = self.geo.phi_foe
        phi1 = phi_s_sa
        ds_SA = self.geo.rb*(0.5*(phi2**2-phi1**2)-self.geo.phi_oo0*(phi2-phi1))
        I = self.CVs.index('sa')
        F[I,:] = ds_SA*self.geo.t/2*(self.p[I,_slice]-p_backpressure)
        
        # Parameters for the S1 chamber
        phi2 = phi_s_sa
        phi1 = phi_s_sa-theta
        ds_S1 = self.geo.rb*(0.5*(phi2**2-phi1**2)-self.geo.phi_oo0*(phi2-phi1))
        I = self.CVs.index('s1')
        F[I,:] = ds_S1*self.geo.t/2*(self.p[I,_slice]-p_backpressure)
        
        # Parameters for the S2 chamber
        phi2 = self.geo.phi_oie
        phi1 = self.geo.phi_oie-theta
        ds_S2 = self.geo.rb*(0.5*(phi2**2-phi1**2)-self.geo.phi_oi0*(phi2-phi1))
        I = self.CVs.index('s2')
        F[I,:] = ds_S2*self.geo.t/2*(self.p[I,_slice]-p_backpressure)
        
        for I in range(1, scroll_geo.nC_Max(self.geo)+1):
            phi2 = self.geo.phi_ooe-pi-theta-2*pi*(I-1)
            phi1 = self.geo.phi_ooe-pi-theta-2*pi*(I)
            ds_C1 = self.geo.rb*(0.5*(phi2**2-phi1**2)-self.geo.phi_oo0*(phi2-phi1))
            ICV = self.CVs.index('c1.'+str(I))
            F[ICV,:] = ds_C1*self.geo.t/2*(self.p[ICV, _slice]-p_backpressure)
            
            phi2 = self.geo.phi_oie-theta-2*pi*(I-1)
            phi1 = self.geo.phi_oie-theta-2*pi*(I)
            ds_C2 = self.geo.rb*(0.5*(phi2**2-phi1**2)-self.geo.phi_oi0*(phi2-phi1))
            ICV = self.CVs.index('c2.'+str(I))
            F[ICV,:] = ds_C2*self.geo.t/2*(self.p[ICV, _slice]-p_backpressure)
        
        phi2 = self.geo.phi_ooe-pi-theta-2*pi*(nC)
        phi1 = self.geo.phi_oos
        ds_D1 = self.geo.rb*(0.5*(phi2**2-phi1**2)-self.geo.phi_oo0*(phi2-phi1))
        ICV = self.CVs.index('d1')
        F[ICV,:] = ds_D1*self.geo.t/2*(self.p[ICV,_slice]-p_backpressure)
        
        phi2 = self.geo.phi_oie-theta-2*pi*(nC)
        phi1 = self.geo.phi_ois
        ds_D2 = self.geo.rb*(0.5*(phi2**2-phi1**2)-self.geo.phi_oi0*(phi2-phi1))
        ICV = self.CVs.index('d2')
        F[ICV,:] = ds_D2*self.geo.t/2*(self.p[ICV,_slice]-p_backpressure)
        
        phi2 = self.geo.phi_ois
        phi1 = self.geo.phi_oi0
        ds_DD = self.geo.rb*(0.5*(phi2**2-phi1**2)-self.geo.phi_oi0*(phi2-phi1))
        ICV = self.CVs.index('dd')
        F[ICV,:] = ds_DD*self.geo.t/2*(self.p[ICV,_slice]-p_backpressure)
        
        ICV = self.CVs.index('ddd')
        F[ICV,:] = (ds_D1+ds_D2+ds_DD)*self.geo.t/2*(self.p[ICV,_slice]-p_backpressure)
        
        # Remove all the nan placeholders
        F[np.isnan(F)] = 0
        
        return F
    
    def attach_HDF5_annotations(self, fName):
        """
        In this function, annotations can be attached to each HDF5 field
        
        Here we add other scroll-specific terms
        
        Parameters
        ----------
        fName : string
            The file name for the HDF5 file that is to be used
        """ 
        #Use the base annotations
        PDSimCore.attach_HDF5_annotations(self, fName)
        
        attrs_dict = {
                '/forces/F_B':'The normal force applied to the orbiting scroll journal bearing [kN]',
                '/forces/F_B0':'The normal force applied to the orbiting scroll journal bearing [kN]',
                '/forces/Fkey':'The forces applied at each key of the Oldham ring [kN]',
                '/forces/Fbackpressure':'The force applied to the orbiting scroll due to the back pressure [kN]',
                '/forces/Fm':'The normal force applied to the orbiting scroll journal bearing [kN]',
                '/forces/Fr':'The radial gas force on the orbiting scroll [kN]',
                '/forces/Ft':'The tangential gas force on the orbiting scroll [kN]',
                '/forces/Fx':'The gas force on the orbiting scroll in the x-direction [kN]',
                '/forces/Fy':'The gas force on the orbiting scroll in the y-direction [kN]',
                '/forces/Fz':'The gas force on the orbiting scroll in the negative z-direction [kN]',
                '/forces/M_B':'The journal bearing moment on the orbiting scroll in the positive x-direction [kN-m]',
                '/forces/M_B0':'The journal bearing moment on the orbiting scroll in the positive x-direction [kN-m]',
                '/forces/Moverturn':'The magnitude of the overturning moment on the orbiting scroll [kN-m]',
                '/forces/Mx':'The overturning moment around the x-axis [kN-m]',
                '/forces/My':'The overturning moment around the y-axis [kN-m]',
                '/forces/Mz':'The spinning moment from the gas around the z-axis [kN-m]',
                '/forces/THETA':'The shifted angle that is used to locate the center of the orbiting scroll [rad]',
                '/forces/Wdot_F1':'The frictional dissipation at key 1 of Oldham ring [kW]',
                '/forces/Wdot_F2':'The frictional dissipation at key 2 of Oldham ring [kW]',
                '/forces/Wdot_F3':'The frictional dissipation at key 3 of Oldham ring [kW]',
                '/forces/Wdot_F4':'The frictional dissipation at key 4 of Oldham ring [kW]',
                '/forces/Wdot_OS_journal':'The frictional dissipation at journal bearing of orbiting scroll [kW]',
                '/forces/Wdot_upper_journal':'The frictional dissipation at upper journal bearing [kW]',
                '/forces/Wdot_lower_journal':'The frictional dissipation at lower journal bearing [kW]',
                '/forces/Wdot_thrust':'The frictional dissipation from thrust bearing [kW]',
                '/forces/Wdot_total':'The frictional dissipation of bearings and friction [kW]',
                '/forces/Wdot_total_mean':'The mean frictional dissipation of bearings and friction over one rotation[kW]',
                '/forces/cx':'The x-coordinate of the centroid of each control volume [m]',
                '/forces/cx_thrust':'Effective x-coordinate of the centroid of all chambers [m]',
                '/forces/cy':'The y-coordinate of the centroid of each control volume [m]',
                '/forces/cy_thrust':'Effective y-coordinate of the centroid of all chambers [m]',
                '/forces/fxp':'Fx/p from geometric analysis for each control volume [kN/kPa]',
                '/forces/fyp':'Fy/p from geometric analysis for each control volume [kN/kPa]',
                '/forces/inertial':'Magnitude of inertial force (m*omega^2*r) [kN]',
                '/forces/mean_Fm':'Mean of Fm over one rotation [kN]',
                '/forces/mean_Fr':'Mean of Fr over one rotation [kN]',
                '/forces/mean_Ft':'Mean of Ft over one rotation [kN]',
                '/forces/mean_Fz':'Mean of Fz over one rotation [kN]',
                '/forces/mean_Mz':'Mean of Mz over one rotation [kN-m]',
                '/forces/mean_tau':'Mean of tau over one rotation [kN-m]',
                '/forces/summed_Fr':'Summation of CV contributions to Fr [kN]',
                '/forces/summed_Ft':'Summation of CV contributions to Ft [kN]',
                '/forces/summed_Fx':'Summation of CV contributions to Fx [kN]',
                '/forces/summed_Fy':'Summation of CV contributions to Fy [kN]',
                '/forces/summed_Fz':'Summation of CV contributions to Fz [kN]',
                '/forces/summed_Moverturn':'Summation of CV contributions to overturning moment [kN-m]',
                '/forces/summed_Mx':'Summation of CV contributions to Mx [kN-m]',
                '/forces/summed_My':'Summation of CV contributions to My [kN-m]',
                '/forces/summed_Mz':'Summation of CV contributions to Mz [kN-m]',
                '/forces/summed_gas_Fz':'Summation of CV contributions to Fz (only from the CV) [kN-m]',
                '/forces/tau':'Torque generated by gas around the central axis of shaft [kN-m]',
                '/forces/xpin':'x-coordinate of the orbiting scroll center [m]',
                '/forces/ypin':'y-coordinate of the orbiting scroll center [m]',
                '/mech/D_crank_bearing':'Diameter of the crank journal bearing [m]',
                '/mech/D_lower_bearing':'Diameter of the lower journal bearing [m]',
                '/mech/D_upper_bearing':'Diameter of the upper journal bearing [m]',
                '/mech/L_crank_bearing':'Length of the crank journal bearing [m]',
                '/mech/L_lower_bearing':'Length of the lowe journal bearing [m]',
                '/mech/L_ratio_bearings':'Ratio of the distances from the upper bearing to the crank bearing [m]',
                '/mech/L_upper_bearing':'Length of the upper journal bearing [m]',
                '/mech/c_crank_bearing':'Clearance (D/2-rb) of the crank journal bearing [m]',
                '/mech/c_lower_bearing':'Clearance (D/2-rb) of the lower journal bearing [m]',
                '/mech/c_upper_bearing':'Clearance (D/2-rb) of the upper journal bearing [m]',
                '/mech/detailed_analysis':'True if detailed mechanical analysis is being used',
                '/mech/journal_tune_factor':'Tuning factor that muliplies losses in each journal bearing [-]',
                '/mech/mu_oil':'Viscosity of the oil [Pa-s]',
                '/mech/oldham_key_friction_coefficient':'Friction coefficient for all keys in Oldham ring [-]',
                '/mech/oldham_key_height':'Height of each key of Oldham ring [m]',
                '/mech/oldham_key_width':'Width of each key of Oldham ring [m]',
                '/mech/oldham_mass':'Mass of Oldham ring [kg]',
                '/mech/oldham_ring_radius':'Radius of Oldham ring [m]',
                '/mech/oldham_rotation_beta':'Angle between Oldham sliding axis and x-axis [radian]',
                '/mech/oldham_thickness':'Thickness of the main ring of Oldham ring [m]',
                '/mech/orbiting_scroll_mass':'Mass of orbiting scroll [kg]',
                '/mech/scroll_density':'Scroll density [kg]',
                '/mech/scroll_plate_thickness':'Scroll plate thickness [m]',
                '/mech/scroll_added_mass':'Scroll added mass [kg]',
                '/mech/scroll_plate_diameter':'Scroll plate diameter [m]',
                '/mech/thrust_ID':'Thrust bearing inner diameter [m]',
                '/mech/thrust_OD':'Thrust bearing outer diameter [m]',
                '/mech/thrust_friction_coefficient':'Thrust bearing friction coefficient [m]'
                }
        
        import h5py
        hf = h5py.File(fName,'a')
        
        for k, v in attrs_dict.iteritems():
            dataset = hf.get(k)
            if dataset is not None:
                dataset.attrs['note'] = v
        hf.close()
        
        
