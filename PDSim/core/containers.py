import copy
import collections
from cPickle import dumps,loads
from CoolProp.State import State

from _containers import TubeCollection, CVArrays
from PDSim.misc.datatypes import arraym
    
class Tube(object):
    """
    A tube is a component of the model that allows for heat transfer and pressure drop.
    
    With this class, the state of at least one of the points is fixed.  For instance, at the inlet of the compressor, the state well upstream is quasi-steady.
    """
    def __init__(self,key1,key2,L,ID,State1=None,State2=None,OD=-1,fixed=-1,TubeFcn=None,mdot=-1,exists=True):
        self.key1 = key1
        self.key2 = key2
        self.fixed = fixed
        
        #: Additional heat to be added to the tube
        self.Q_add = 0.0
        
        #: Fixed heat transfer coefficient if desired (if less than zero will use correlation - default)
        self.alpha = -1.0
        
        self.exists = exists
        if fixed<0:
            raise AttributeError(textwrap.dedent("""You must provide an integer 
            value for fixed, either 1 for Node 1 fixed, or 2 for Node 2 fixed.  
            You provided None (or didn\'t include the parameter"""))
        if fixed==1 and isinstance(State1,State) and State2==None:
            #Everything good
            self.State1=State1
            self.State2=State(self.State1.Fluid,{'T':self.State1.T,'D':self.State1.rho})
        elif fixed==2 and isinstance(State2,State) and State1==None:
            #Everything good
            self.State2=State2
            self.State1=State(self.State2.Fluid,{'T':self.State2.T,'D':self.State2.rho})
        else:
            raise AttributeError('Incompatibility between the value for fixed and the states provided')
            
        self.TubeFcn=TubeFcn
        if mdot<0:
            self.mdot=0.010
            print('Warning: mdot not provided to Tube class constructor, guess value of '+str(self.mdot)+' kg/s used')
        else:
            self.mdot=mdot
        self.L=L
        self.ID=ID
        self.OD=OD

from _containers import ControlVolume, ControlVolumeCollection
