
class PDSimPlugin(object):
    """
    This is the base class that represents plugins for the GUI
    """
     
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
        

    def _check_plugin_chunks(self, chunks):
        """ Check whether the script chunks are valid """
        
        allowed = ['pre_import','post_import','pre_build','pre_build_instantiation','post_build_instantiation','post_build','pre_run','post_run','plugin_injected_chunks']
        for key in chunks.keys():
            if key not in allowed:
                raise ValueError(key+' is an invalid key to be returned in dictionary from get_script_chunks')
        
    def get_script_chunks(self):
        """
        Get the chunks for the script from the plugin
        
        Return a dictionary of strings for the script, the keys that are allowed are:
        
        * ``pre_import`` (goes before all the standard imports)
        * ``post_import`` (goes after all the standard imports)
                
        * ``pre_build`` (goes at the very beginning of the build function)
        * ``pre_build_instantiation`` (goes right before instantiation of the PDSimCore subclass)
        * ``post_build_instantiation`` (goes right before instantiation of the PDSimCore subclass)
        * ``post_build`` (goes at the very end of the build function)
        
        
        * ``pre_run`` (goes at the beginning of the run function)
        * ``post_run`` (goes at the end of the run function)
        
        * ``plugin_injected_chunks`` (a dictionary of chunks that is passed to the ``InputToolBook.get_script_chunks()`` function 
        
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
    
    
