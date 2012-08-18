
class PDSimPlugin(object):
    def __init__(self):
        """
        A base class that is to be subclassed for other plugins to PDSim.
        
        Plugins are structured in order to be as flexible as possible.
        """
        pass
    
    def apply(self, sim, **kwargs):
        """
        Apply the plugin's code - it can do whatever it wants to the simulation 
        """
        raise NotImplementedError("Subclasses of PDSimPlugin must provide the apply() function")
        
    def set_GUI(self, Main):
        self.GUI = Main
        
    def activate(self, event):
        self._activated = True
        
    def deactivate(self):
        self._activated = True
        
    def is_activated(self):
        return self._activated
    
