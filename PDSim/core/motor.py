import operator

# If scipy is available, use its spline interpolation function, otherwise, 
# use our implementation (for packaging purposes)
try:
    import scipy.interpolate as interp
except ImportError:
    import PDSim.misc.scipylike as interp
    
import warnings

class Motor(object):
    """
    A class that implements the motor model
    
    This can be either a constant efficiency model
    """
    
    def __init__(self):
        
        self.type = ''
        self.suction_fraction = 0.0
            
    def set_eta(self, eta):
        """
        Set the motor efficiency
        """
        self.eta_motor = eta
        self.type = 'const_eta_motor'
        
    def set_coeffs(self, tau_coeffs, eta_coeffs, omega_coeffs):
        """
        Set the coefficients
        
        Parameters
        ----------
        tau_coeffs: iterable (list or similar)
            Values for the torque [N-m]
        eta_coeffs: iterable (list or similar)
            Values for the efficiency [-] (in the range 0 to 1)
        omega_coeffs: iterable (list or similar)
            Values for the rotational speed [rad/s]
        """
        #Only accept coefficients if they are all the same length
        assert len(tau_coeffs) == len(eta_coeffs) == len(omega_coeffs)
        #Set the local values
        self.tau_coeffs = tau_coeffs
        self.eta_coeffs = eta_coeffs
        self.omega_coeffs = omega_coeffs
        
        #interp1d needs the independent variable to be monotonically increasing
        #Join the lists together 
        zipped = zip(self.tau_coeffs, self.eta_coeffs, self.omega_coeffs)
        #Sort the joined lists by the torque
        zipped = sorted(zipped, key = operator.itemgetter(0))
        #Unzip them back into the original variables
        self.tau_coeffs, self.eta_coeffs, self.omega_coeffs = zip(*zipped)
        self.type = 'motor_map'
        
    def plot_eta(self, ax):
        """
        Plot eta v. torque on the given axis
        """
        pass
    
    def plot_speed(self, ax):
        """
        Plot slip speed v. torque on the given axis
        """
        pass
     
    def invert_map(self, Wdot, kind = 'linear'):
        """
        Invert the map to calculate the speed and the torque based on the power
        the power is given by tau*omega
        
        If a constant efficiency, just return (efficiency,None) tuple
        """
        if not kind == 'linear':
            warnings.warn('invert_map does not take parameter "kind" anymore')
            
        if self.type == 'const_eta_motor':
            return self.eta_motor,None
        else:
            Wdot_coeffs = [tau*omega/1000 for tau,omega in zip(self.tau_coeffs, self.omega_coeffs)]
            #Do the 1D interpolation
            eta_interp = interp.splrep(Wdot_coeffs, self.eta_coeffs, k=2, s=0)
            eta = interp.splev(Wdot, eta_interp)
            omega_interp = interp.splrep(Wdot_coeffs, self.omega_coeffs, k=2, s=0)
            omega = interp.splev(Wdot, omega_interp)
        
            return eta,omega
    
    def apply_map(self, tau, kind = 'linear'):
        """        
        Actually use the motor map to calculate the slip speed and the motor
        efficiency
        
        Parameters
        ----------
        tau : float
            Torque [N-m]
        kind : string, optional
            The kind of interpolation to do (see scipy.interpolate.interp1d) - deprecated
            
        Returns
        -------
        eta : float
            Efficiency [-]
        omega : float
            Rotational speed [rad/s]
        """
        if not kind == 'linear':
            warnings.warn('apply_map does not take parameter "kind" anymore')
            
        eta_interp = interp.splrep(self.tau_coeffs, self.eta_coeffs, k=2, s=0)
        eta = interp.splev(tau, eta_interp)
        omega_interp = interp.splrep(self.tau_coeffs, self.omega_coeffs, k=2, s=0)
        omega = interp.splev(tau, omega_interp)
        
        #Return the values
        return eta, omega