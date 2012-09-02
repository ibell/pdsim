
class PDSimPlugin(object):
    def __init__(self):
        """
        A base class that is to be subclassed for other plugins to PDSim.
        
        Plugins are structured in order to be as flexible as possible.
        """
        self._activated = False
        self.GUI = None
        
    def should_enable(self):
        """
        Pre-enabling code that determines whether the plugin should be enabled.
        
        Can be overloaded in the derived class
        """
        return True
    
    def apply(self, sim, **kwargs):
        """
        Apply the plugin's code - it can do whatever it wants to the simulation 
        """
        raise NotImplementedError("Subclasses of PDSimPlugin must provide the apply() function")
        
    def set_GUI(self, Main):
        self.GUI = Main
        
    def activate(self, event):
        """
        Activate the plugin
        """
        self._activated = not self._activated
        
    def is_activated(self):
        return self._activated
    
    def post_process(self, simulation):
        """
        Do any post-processing required of the data after the model has been
        run to completion
        """
        pass

    def collect_output_terms(self):
        return []
    
    
