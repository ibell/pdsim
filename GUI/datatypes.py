from PDSim.misc.datatypes import AnnotatedValue
import wx

def HeaderStaticText(parent, caption):
    """
    Make a static text with black background and white text
    """
    text = wx.StaticText(parent, label = caption)
    text.SetBackgroundColour((0, 0, 0)) # set text color
    text.SetForegroundColour((255, 255, 255)) # set text back color
    return text
    
class AnnotatedGUIObject(AnnotatedValue):
    def __init__(self, obj, GUI_location):
        """
        A wrapper around the AnnotatedValue class
        Parameters
        ----------
        obj : :class:`AnnotatedValue <PDSim.misc.datatypes.AnnotatedValue>` instance
        GUI_location : The wxpython element in the GUI that stores the value of this object.
        """
        for attr in dir(obj):
            # Not a private method
            if attr.find('__')<0:
                # Copy each field over from AnnotatedValue
                setattr(self, attr, getattr(obj,attr))
                
        self.GUI_location = GUI_location
        if not hasattr(self.GUI_location,'GetValue'):
            raise AttributeError('Could not create AnnotatedGUIObject, GUI_location [{t:s}] does not implement GetValue'.format(t=self.GUI_location))
        if not hasattr(self.GUI_location,'SetValue'):
            raise AttributeError('Could not create AnnotatedGUIObject, GUI_location [{t:s}] does not implement SetValue'.format(t=self.GUI_location))
        
    def GetValue(self):
        return self.GUI_location.GetValue()
    
    def SetValue(self, value):
        self.GUI_location.SetValue(value)
        
class CoupledAnnotatedGUIObject(AnnotatedGUIObject):
    """
    This is an annotated gui object that must be set in partnership with 
    another one or more terms. For instance, the inlet state to the machine 
    requires two state parameters.
    
    Example: If you change the saturation suction temperature, it is not clear 
    how you intend to change the superheat, so you must provide both
    
    Parameters
    ----------
    obj : :class:`AnnotatedValue <PDSim.misc.datatypes.AnnotatedValue>` instance
    GUI_location : The wxpython element in the GUI that stores the value of this object.
    handler : A function that is called with the list of linked parameters
    
    """
    
    def __init__(self, 
                 obj, 
                 GUI_location,  
                 handler
                 ):
        AnnotatedGUIObject.__init__(self, obj, GUI_location)
        self.handler = handler
        
    def link_required_parameters(self, required_partners):
        """
        Parameters
        ----------
        required_partners : list, optional
            A list of other :class:`AnnotatedGUIObjects` that must be provided for this object
            The other linked objects must also link back to this object
        """
        if isinstance(required_partners, list):
            self.required_partners = required_partners
        else:
            self.required_partners = [required_partners]
        for partner in self.required_partners:
            if not isinstance(partner, CoupledAnnotatedGUIObject):
                raise KeyError('partner ({o:s}) not a CoupledAnnotatedGUIObject'.format(o = partner))

class InfiniteList(object):
    """
    Creates a special list where removing an element just puts it back at the end of the list
    """
    def __init__(self, values):
        """
        Parameters
        ----------
        values : list

        """
        self.values = values
        
    def pop(self):
        """
        Return the first element, then put the first element back at the end of the list
        """
        val1 = self.values[0]
        self.values.pop(0)
        self.values.append(val1)
        return val1
    
    def prepend(self,item):
        """
        Put the item back to the front of the list
        
        Parameters
        ----------
        item : object
            Thing to push to the beginning of the list
            
        """
        self.values.remove(item)
        self.values.insert(0,item)