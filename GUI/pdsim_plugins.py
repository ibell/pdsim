import abc

class PDSimPlugin(object):
    __metaclass__ = abc.ABCMeta
     
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
    
    @abc.abstractmethod
    def apply(self, **kwargs):
        """
        Apply the plugin's code - it can do whatever it wants to the GUI 
        """
        raise NotImplementedError("Subclasses of PDSimPlugin must provide the apply() function")
        

    def _check_script_chunks(self):
        """ Check whether the script chunks are valid """
        
        allowed = ['post_imports','pre_build','post_build']
        return all([key in allowed for key in chunks.keys()])
    
    def get_script_chunks(self):
        """
        Get the chunks for the script from the plugin
        
        Return a dictionary of chunks (strings) for the script, the keys that are allowed are:
        
        * ``post_imports`` (goes after all the standard imports)
        * ``pre_build`` (goes at the beginning of the build function)
        * ``post_build`` (goes at the end of the build function)
        
        """
        return {}

    def set_GUI(self, Main):
        self.GUI = Main
        
    def activate(self, event):
        """
        Function to activate the plugin
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
    
    
