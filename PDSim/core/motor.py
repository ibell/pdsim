import operator
from scipy.interpolate import interp1d

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
        
        if self.type == 'const_eta_motor':
            return self.eta_motor,None
        else:
            Wdot_coeffs = [tau*omega/1000 for tau,omega in zip(self.tau_coeffs, self.omega_coeffs)]
            #Do the 1D interpolation
            eta = float(interp1d(Wdot_coeffs, self.eta_coeffs, kind = kind)(Wdot))
            omega = float(interp1d(Wdot_coeffs, self.omega_coeffs, kind = kind)(Wdot))
            
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
            The kind of interpolation to do (see scipy.interpolate.interp1d)
            
        Returns
        -------
        eta : float
            Efficiency [-]
        omega : float
            Rotational speed [rad/s]
        """
        
        #Do the 1D interpolation
        eta = float(interp1d(self.tau_coeffs, self.eta_coeffs, kind = kind)(tau))
        omega = float(interp1d(self.tau_coeffs, self.omega_coeffs, kind = kind)(tau))
        
        #Return the values
        return eta, omega
        
        