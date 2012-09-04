# -*- coding: latin-1 -*-

import wx
from wx.lib.mixins.listctrl import CheckListCtrlMixin
import CoolProp
from CoolProp.State import State
from CoolProp import CoolProp as CP
from ConfigParser import SafeConfigParser
import codecs
import numpy as np
import os
import itertools
import matplotlib as mpl
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as WXCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as WXToolbar
from multiprocessing import Process
from PDSim.scroll import scroll_geo

import quantities as pq

length_units = {
                'Meter': pq.length.m,
                'Micrometer' : pq.length.um,
                'Centimeter' : pq.length.cm,
                'Inch' : pq.length.inch,
                }

area_units = {
                'Square Meter': pq.length.m**2,
                'Square Micrometer' : pq.length.um**2,
                'Square Centimeter' : pq.length.cm**2,
                'Square Inch' : pq.length.inch**2,
                }

volume_units = {
                'Cubic Meter': pq.length.m**3,
                'Cubic Micrometer' : pq.length.um**3,
                'Cubic Centimeter' : pq.length.cm**3,
                'Cubic Inch' : pq.length.inch**3,
                }

pressure_units = {
                  'kPa' : pq.kPa,
                  'psia' : pq.psi
                  }
rev = pq.UnitQuantity('revolution', 2*np.pi*pq.radians, symbol='rev')

rotational_speed_units ={
                         'Radians per second': pq.radians/pq.sec,
                         'Radians per minute': pq.radians/pq.min,
                         'Revolutions per second': rev/pq.sec,
                         'Revolutions per minute': rev/pq.min,
                         }

temperature_units = {
                     'Kelvin' : np.nan,
                     'Celsius' : np.nan,
                     'Fahrenheit' : np.nan,
                     'Rankine': np.nan
                     }

class UnitConvertor(wx.Dialog):
    def __init__(self, value, default_units, type = None, TextCtrl = None):
        wx.Dialog.__init__(self, None, title='Convert units')
        
        self.default_units = default_units
        self.__is_temperature__ = False
        if default_units in length_units or type == 'length':
            self.unit_dict = length_units
        elif default_units in area_units or type == 'area':
            self.unit_dict = area_units
        elif default_units in volume_units or type == 'volume':
            self.unit_dict = volume_units
        elif default_units in rotational_speed_units or type == 'rotational_speed':
            self.unit_dict = rotational_speed_units
        elif default_units in pressure_units or type == 'pressure':
            self.unit_dict = pressure_units
        elif default_units in temperature_units or type == 'temperature':
            self.unit_dict = temperature_units
            self.__is_temperature__ = True
        else:
            raise KeyError('Sorry your units '+default_units+' did not match any of the unit terms')
            
        self.txt = wx.TextCtrl(self, value=str(value))
        self.units = wx.ComboBox(self)
        self.units.AppendItems(sorted(self.unit_dict.keys()))
        self.units.SetEditable(False)
        if default_units in self.units.GetStrings():
            self.units.SetStringSelection(default_units)
        else:
            raise KeyError
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.txt)
        sizer.Add(self.units)
        self.SetSizer(sizer)
        self.Fit()
        self._old_units = default_units
        
        self.Bind(wx.EVT_COMBOBOX, self.OnSwitchUnits, self.units)
        
    def OnSwitchUnits(self, event):
        if not self.__is_temperature__:
            old = float(self.txt.GetValue()) * self.unit_dict[self._old_units]
            new_units_str = self.units.GetStringSelection()
            new_units = self.unit_dict[new_units_str]
            old.units = new_units
            self._old_units = new_units_str
            self.txt.SetValue(str(old.magnitude))
        else:
            old_val = float(self.txt.GetValue())
            new_units_str = self.units.GetStringSelection()
            new_val = self._temperature_convert(self._old_units, old_val, new_units_str)
            self._old_units = new_units_str
            self.txt.SetValue(str(new_val))
            
    def get_value(self):
        """
        Return a string with value in the original units for use in calling function to dialog
        """
        if not self.__is_temperature__:
            old = float(self.txt.GetValue()) * self.unit_dict[self._old_units]
            old.units = self.unit_dict[self.default_units]
            return str(old.magnitude)
        else:
            old_val = float(self.txt.GetValue())
            return str(self._temperature_convert(self._old_units, old_val, self.default_units))
    
    def _temperature_convert(self, old, old_val, new):
        """
        Internal method to convert temperature
        
        Parameters
        ----------
        old : string
        old_val : float
        new : string
        """
        #convert old value to Celsius
        # also see: http://en.wikipedia.org/wiki/Temperature
        if old == 'Fahrenheit':
            celsius_val = (old_val-32)*5.0/9.0
        elif old == 'Kelvin':
            celsius_val = old_val-273.15
        elif old == 'Rankine':
            celsius_val = (old_val-491.67)*5.0/9.0
        elif old == 'Celsius':
            celsius_val = old_val
            
        #convert celsius to new value
        if new == 'Celsius':
            return celsius_val
        elif new == 'Fahrenheit':
            return celsius_val*9.0/5.0+32.0
        elif new == 'Kelvin':
            return celsius_val+273.15
        elif new == 'Rankine':
            return (celsius_val+273.15)*9.0/5.0
        
        
class PlotPanel(wx.Panel):
    def __init__(self, parent, **kwargs):
        wx.Panel.__init__(self, parent, size = (300,200), **kwargs)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.figure = mpl.figure.Figure(dpi=100, figsize=(2, 2))
#        self.figure.set_figwidth(2.0)
#        self.figure.set_figheight(2.0)
        self.canvas = WXCanvas(self, -1, self.figure)
#        self.canvas.resize(200,200)
        self.toolbar = WXToolbar(self.canvas)
        self.toolbar.Realize()
        sizer.Add(self.canvas)
        sizer.Add(self.toolbar)
        self.SetSizer(sizer)
        sizer.Layout()

class PDPanel(wx.Panel):
    """
    A base class for panel with some goodies thrown in 
    
    Not intended for direct instantiation, rather it should be 
    subclassed, like in :class:`recip_panels.GeometryPanel`
    
    Loading from configuration file
    -------------------------------
    Method (A)
    
    Any subclassed PDPanel must either store a list of items as 
    ``self.items`` where the entries are dictionaries with at least the entries:
    
    - ``text``: The text description of the term 
    - ``attr``: The attribute in the simulation class that is linked with this term
    
    Other terms that can be included in the item are:
    - ``tooltip``: The tooltip to be attached to the textbox
    
    Method (B)
    
    The subclassed panel can provide the functions post_prep_for_configfile
    and post_get_from_configfile, each of which take no inputs.
    
    In post_prep_for_configfile, the subclassed panel can package up its elements
    into a form that can then be re-created when post_get_from_configfile is
    called
    
    
    Saving to configuration file
    ----------------------------
    
    Adding terms to parametric table
    --------------------------------
    
    """
    def __init__(self,*args,**kwargs):
        wx.Panel.__init__(self,*args,**kwargs)
        self.name=kwargs.get('name','')
        
    def _get_item_by_attr(self, attr):
        if hasattr(self,'items'):
            for item in self.items:
                if item['attr'] == attr:
                    return item
        raise ValueError('_get_item_by_attr failed')
        
    def _get_value(self,thing):
        #This first should work for wx.TextCtrl
        if hasattr(thing,'GetValue'):
            value=str(thing.GetValue()).strip()
            try:
                return float(value)
            except ValueError:
                return value
        elif hasattr(thing,'GetSelectionString'):
            value=thing.GetSelection()
            try:
                return float(value)
            except ValueError:
                return value
             
    def set_params(self, sim):
        
        if not hasattr(self,'items'):
            return
        else:
            items = self.items
        
        if hasattr(self,'skip_list'):
            # Don't actually set these attributes (they might over-write 
            # methods or attributes in the simulation)
            items = [item for item in items if item['attr'] not in self.skip_list()]
            
        for item in items:
            setattr(sim, item['attr'],self._get_value(item['textbox']))
    
    def ConstructItems(self,items,sizer,configdict=None, descdict=None, parent = None):
        for item in items:
            #item is a dictionary of values including the keys:
            #  - attr
            #  - textbox
            #  - val
            
            if parent is None:
                parent = self
            
            if 'val' not in item and configdict is not None:
                k = item['attr']
                if k not in configdict:
                    self.warn_unmatched_attr(k)
                    val,item['text']=self.get_from_configfile(self.name, k, default = True)
                    print val, item['text']
                else:
                    val = configdict[k]
                    item['text'] = descdict[k]
            else:
                val = item['val']
                item['text'] = descdict[k]
                
            label=wx.StaticText(parent, -1, item['text'])
            sizer.Add(label, 1, wx.EXPAND)
            textbox=wx.TextCtrl(parent,-1,str(val))
            sizer.Add(textbox, 1, wx.EXPAND)
            item.update(dict(textbox=textbox,label=label))
            
            caption = item['text']
            if caption.find(']')>=0 and caption.find(']')>=0: 
                units = caption.split('[',1)[1].split(']',1)[0]
                if units == 'm':
                    textbox.default_units = 'Meter'
                elif units == 'm²':
                    textbox.default_units = 'Square Meter'
                elif units == 'm³':
                    textbox.default_units = 'Cubic Meter'
                elif units == 'rad/s':
                    textbox.default_units = 'Radians per second'
                self.Bind(wx.EVT_CONTEXT_MENU,self.OnChangeUnits,textbox)      
        
    def warn_unmatched_attr(self, attr):
        print "didn't match attribute", attr
        
    def prep_for_configfile(self):
        """
        Writes the panel to a format ready for writing to config file
        using the entries in ``self.items``.  
        
        If there are other fields that need to get saved to file, the panel 
        can provide the ``post_prep_for_configfile`` function and add the additional fields 
        
        This function will call the ``post_prep_for_configfile`` if the subclass has it
        and add to the returned string
        """
        if self.name=='':
            return ''
            
        if not hasattr(self,'items'):
            self.items=[]
        
        s='['+self.name+']\n'
        
        for item in self.items:
            val = item['textbox'].GetValue()
            # Description goes into the StaticText control
            try:
                int(val)
                type_='int'
            except ValueError:
                try: 
                    float(val)
                    type_='float'
                except ValueError:
                    type_='string'
            s+=item['attr']+' = '+type_+','+item['textbox'].GetValue().encode('latin-1')+','+item['text']+'\n'
            
        if hasattr(self,'post_prep_for_configfile'):
            s+=self.post_prep_for_configfile()
        
        s=s.replace('%','%%')
        return s
           
    def _get_from_configfile(self, name, value):
        
        #Split at the first comma to get type, and value+description
        type,val_desc = value.split(',',1)
        #If it has a description, use it, otherwise, just use the config file key
        if len(val_desc.split(','))==2:
            val,desc_=val_desc.split(',')
            desc=desc_.strip()
        else:
            val=val_desc
            desc=name.strip()
            
        if type=='int':
            d=int(val)
        elif type=='float':
            d=float(val)
        elif type=='str':
            d=unicode(val)
        elif type=='State':
            Fluid,T,rho=val.split(',')
            d=dict(Fluid=Fluid,T=float(T),rho=float(rho))
        else:
            #Try to let the panel use the (name, value) directly
            if hasattr(self,'post_get_from_configfile'):
                d = self.post_get_from_configfile(name, value)
            else:
                raise KeyError('Type in line '+name+' = ' +value+' must be one of int,float,str')     
        return d, desc 
               
    def get_from_configfile(self, section, key = None, default = False):
        """
        configfile: file path or readable object (StringIO instance or file)
        Returns a dictionary with each of the elements from the given section 
        name from the given configuration file.  Each of the values in the configuration 
        file may have a string in the format 
        
        int,3,
        float,3.0
        string,hello
        
        so that the code can know what type the value is.  If the value is something else,
        ask post_get_from_configfile if it knows what to do with it
        
        """
        d, desc={}, {}
        
        Main = wx.GetTopLevelParent(self)
        parser, default_parser = Main.get_config_objects()
        
        #Section not included, use the default section from the default config file
        if not parser.has_section(section):
            dlg = wx.MessageDialog(None,'Section '+section+' was not found, falling back to default configuration file')
            dlg.ShowModal()
            dlg.Destroy()
            _parser = default_parser
        elif default:
            # We are programmatically using the default parameters, 
            # don't warn automatically'
            _parser = default_parser
        else:
            _parser = parser
        
        if key is not None and default:
            value = _parser.get(section, key)
            _d,_desc =  self._get_from_configfile(key, value)
            return _d,_desc
        
        for name, value in _parser.items(section):
            _d,_desc = self._get_from_configfile(name,value)
            d[name] = _d
            desc[name] = _desc
            
        return d,desc
    
    def get_additional_parametric_terms(self):
        """
        Provide additional parametric terms to the parametric table builder
        
        This function, when implemented in the derived class, will provide 
        additional terms to the parametric table builder.  If unimplemented, will
        just return ``None``.
        
        The format of the returned list should be a list of dictionaries with
        the terms:
            
            parent : a reference to this panel
            attr : the attribute in the simulation that will be linked with this term
            text : the label for the term
            
        Notes
        -----
        ``attr`` need not be actually the term that will ultimately be set in the simulation model
        
        It will be parsed by the apply_additional_parametric_term() function 
        below
        """
        pass
    
    def apply_additional_parametric_terms(self, attrs, vals, items):
        """
        
        Returns
        -------
        list of items that were unmatched
        
        Raises
        ------
        """
        return attrs, vals
    
    def BindChangeUnits(self, TextCtrl):
        self.Bind(wx.EVT_KEY_DOWN, self.OnChangeUnits, TextCtrl)
        
    def OnChangeUnits(self, event):
        TextCtrl = event.GetEventObject()
        dlg = UnitConvertor(value = float(TextCtrl.GetValue()),
                            default_units = TextCtrl.default_units
                            )
        
        dlg.ShowModal()
        TextCtrl.SetValue(dlg.get_value())
        dlg.Destroy()
        
class ChangeParamsDialog(wx.Dialog):
    def __init__(self, params, **kwargs):
        wx.Dialog.__init__(self, None, **kwargs)
    
        self.params = params
        sizer = wx.FlexGridSizer(cols = 2)
        self.labels = []
        self.values = []
        self.attrs = []
    
        for p in self.params:
            l, v = LabeledItem(self,
                               label = p['desc'],
                               value = str(p['value'])
                               )
            self.labels.append(l)
            self.values.append(v)
            self.attrs.append(p['attr'])
            
            sizer.AddMany([l,v])
            
        self.SetSizer(sizer)
        min_width = min([l.GetSize()[0] for l in self.labels])
        for l in self.labels:
            l.SetMinSize((min_width,-1))
        sizer.Layout()
        self.Fit()
        
         #Bind a key-press event to all objects to get Esc 
        children = self.GetChildren()
        for child in children:
            child.Bind(wx.EVT_KEY_UP,  self.OnKeyPress)
    
    def get_values(self):
        params = []
        for l,v,k in zip(self.labels, self.values, self.attrs):
            params += [dict(desc = l.GetLabel(),
                           attr = k,
                           value = float(v.GetValue())
                        )]
        return params
            
    def OnAccept(self, event = None):
        self.EndModal(wx.ID_OK)
        
    def OnKeyPress(self,event = None):
        """ cancel if Escape key is pressed or accept if Enter """
        event.Skip()
        if event.GetKeyCode() == wx.WXK_ESCAPE:
            self.EndModal(wx.ID_CANCEL)
        elif event.GetKeyCode() == wx.WXK_RETURN:
            self.EndModal(wx.ID_OK)
    
    def CancelValues(self, event = None):
        self.EndModal(wx.ID_CANCEL)
            
class MassFlowOption(wx.Panel):
    def __init__(self,
                 parent,
                 key1, 
                 key2,
                 label = 'NONE',
                 types = None,
                 ):
        """
        A wx.Panel for selecting the flow model and associated parameters
        
        Should not be instantiated directly, rather subclassed in order to provide the list of dictionaries
        of flow models for a given type of machine
        """
        wx.Panel.__init__(self, parent)
        
        self.key1 = key1
        self.key2 = key2
        
        options = self.model_options()
        
        self.label = wx.StaticText(self, label=label)
        self.choices = wx.ComboBox(self)
        
        for option in options:
            self.choices.Append(option['desc'])
        self.choices.SetSelection(0)
        self.choices.SetEditable(False)
        self.options_list = options
        
        self.params = wx.Button(self, label='Params...')
        if not 'params' in option or not option['params']:
            self.params.Enable(False)
            
        else:
            TTS = self.dict_to_tooltip_string(option['params'])
            self.params.SetToolTipString(TTS)
            self.params_dict = option['params']
            
            self.params.Bind(wx.EVT_BUTTON, self.OnChangeParams)
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.label)
        sizer.Add(self.choices)
        sizer.Add(self.params)
        self.SetSizer(sizer)
        sizer.Layout()
    
    def OnChangeParams(self, event):
        """
        Open a dialog to change the values
        """
        dlg = ChangeParamsDialog(self.params_dict)
        if dlg.ShowModal() == wx.ID_OK:
            self.params_dict = dlg.get_values()
            TTS = self.dict_to_tooltip_string(self.params_dict)
            self.params.SetToolTipString(TTS)
        dlg.Destroy()
        
    def dict_to_tooltip_string(self, params):
        s = ''
        for param in params:
            s += param['desc'] + ': ' + str(param['value']) + '\n'
        return s
    
    def get_function_name(self):
        for option in self.options_list:
            if option['desc'] == self.choices.GetStringSelection():
                return option['function_name']
        raise AttributeError
        
    def model_options(self):
        """
        This function should return a list of dictionaries.  
        In each dictionary, the following terms must be defined:
        
        * desc : string
            Very terse description of the term 
        * function_name : function
            the function in the main machine class to be called
        * params : list of dictionaries
        
        MUST be implemented in the sub-class
        """
        raise NotImplementedError
    
class ParaSelectDialog(wx.Dialog):
    def __init__(self):
        wx.Dialog.__init__(self, None, title = "State Chooser",)
        self.MinLabel, self.Min = LabeledItem(self, label = 'Minimum value', value = '')
        self.MaxLabel, self.Max = LabeledItem(self, label = 'Maximum value', value = '')
        self.StepsLabel, self.Steps = LabeledItem(self, label = 'Number of steps', value = '')
        self.Accept = wx.Button(self, label = "Accept")
        sizer = wx.FlexGridSizer(cols = 2, hgap = 4, vgap = 4)
        sizer.AddMany([self.MinLabel, self.Min, self.MaxLabel, 
                       self.Max, self.StepsLabel, self.Steps, self.Accept])
        self.SetSizer(sizer)
        sizer.Layout()
        self.Fit()
        self.Accept.Bind(wx.EVT_BUTTON, self.OnAccept)
        
        #Bind a key-press event to all objects to get Esc 
        children = self.GetChildren()
        for child in children:
            child.Bind(wx.EVT_KEY_UP,  self.OnKeyPress)
        
    def join_values(self):
        values = np.linspace(float(self.Min.GetValue()),float(self.Max.GetValue()),int(self.Steps.GetValue()))
        return ', '.join([str(val) for val in values]) 
        
    def OnAccept(self, event = None):
        self.EndModal(wx.ID_OK)
        
    def OnKeyPress(self,event = None):
        """ cancel if Escape key is pressed """
        event.Skip()
        if event.GetKeyCode() == wx.WXK_ESCAPE:
            self.EndModal(wx.ID_CANCEL)
    
    def CancelValues(self, event = None):
        self.EndModal(wx.ID_CANCEL)

class ParametricOption(wx.Panel):
    def __init__(self, parent, items):
        wx.Panel.__init__(self, parent)
        
        attrs = [item['attr'] for item in items]
        labels = [item['text'] for item in items]
        self.Terms = wx.ComboBox(self)
        self.Terms.AppendItems(labels)
        self.Terms.SetSelection(0)
        self.Terms.SetEditable(False)
        self.RemoveButton = wx.Button(self, label = '-', style = wx.ID_REMOVE)
        self.Values = wx.TextCtrl(self, value = '1,2,3,4,5,6,7,8,9')
        self.Select = wx.Button(self, label = 'Select...')
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.RemoveButton)
        sizer.Add(self.Terms)
        sizer.Add(self.Values)
        sizer.Add(self.Select)
        self.SetSizer(sizer)
        self.Select.Bind(wx.EVT_BUTTON,self.OnSelectValues)
        self.RemoveButton.Bind(wx.EVT_BUTTON, lambda event: self.Parent.RemoveTerm(self))
    
    def OnSelectValues(self, event = None):
        dlg = ParaSelectDialog()
        if dlg.ShowModal() == wx.ID_OK:
            self.Values.SetValue(dlg.join_values())
        dlg.Destroy()
        
    def get_values(self):
        name = self.Terms.GetStringSelection()
        #To list of floats
        values = [float(val) for val in self.Values.GetValue().split(',')]
        return name, values
    
    def set_values(self,key,value):
        self.Terms.SetStringSelection(key)
        self.Values.SetValue(value)
        
    def update_parametric_terms(self, items):
        """
        Update the items in each of the comboboxes
        """
        labels = [item['text'] for item in items]
        #Get the old string
        old_val = self.Terms.GetStringSelection()
        #Update the contents of the combobox
        self.Terms.SetItems(labels)
        #Reset the string
        self.Terms.SetStringSelection(old_val)
        
        
class ParametricCheckList(wx.ListCtrl, CheckListCtrlMixin):
    def __init__(self, parent, headers, values):
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        CheckListCtrlMixin.__init__(self)
        
        #Build the headers
        self.InsertColumn(0, '')
        self.SetColumnWidth(0, wx.LIST_AUTOSIZE)
        for i, header in enumerate(headers):
            self.InsertColumn(i+1, header)
        
        self.data = [row for row in itertools.product(*values)]
        
        #Add the values one row at a time
        for i,row in enumerate(self.data):
            self.InsertStringItem(i,'')
            for j,val in enumerate(row):
                self.SetStringItem(i,j+1,str(val))
            self.CheckItem(i)
            
        for i in range(len(headers)):
            self.SetColumnWidth(i+1,wx.LIST_AUTOSIZE_USEHEADER)
        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)

    def OnItemActivated(self, event):
        self.ToggleItem(event.m_itemIndex)
    
    def GetStringItem(self,Irow,Icol):
        return self.data[Irow][Icol]
                                 
class ParametricPanel(PDPanel):
    def __init__(self, parent, configfile, items, **kwargs):
        PDPanel.__init__(self, parent, **kwargs)
        
        self.variables =  items
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.ButtonSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.AddButton = wx.Button(self, label = "Add Term", style = wx.ID_ADD)
        self.AddButton.Bind(wx.EVT_BUTTON, self.OnAddTerm)
        self.ButtonSizer.Add(self.AddButton)
        sizer.Add(self.ButtonSizer)
        self.SetSizer(sizer)
        sizer.Layout()
        self.NTerms = 0
        self.ParamSizer = None
        self.ParamListSizer = None
        self.ParaList = None
        self.RunButton = None
        
        #Has no self.items, so all processing done through post_get_from_configfile
        self.get_from_configfile('ParametricPanel')
        
    def OnAddTerm(self, event = None):
        if self.NTerms == 0:
            self.ParamSizerBox = wx.StaticBox(self, label = "Parametric Terms")
            self.ParamSizer = wx.StaticBoxSizer(self.ParamSizerBox, wx.VERTICAL)
            self.GetSizer().Add(self.ParamSizer)
            self.BuildButton = wx.Button(self, label = "Build Table")
            self.BuildButton.Bind(wx.EVT_BUTTON, self.OnBuildTable)
            self.ButtonSizer.Add(self.BuildButton)
        option = ParametricOption(self, self.variables)
        self.ParamSizer.Add(option)
        self.ParamSizer.Layout()
        self.NTerms += 1
        self.GetSizer().Layout()
        self.Refresh()
    
    def RemoveTerm(self, term):
        term.Destroy()
        self.NTerms -= 1
        if self.NTerms == 0:
            self.BuildButton.Destroy()
            if self.ParamSizer is not None:
                self.GetSizer().Remove(self.ParamSizer)
                self.ParamSizer = None
            if self.ParaList is not None:
                self.ParaList.Destroy()
                self.ParaList = None
            if self.ParamListSizer is not None:
                self.GetSizer().Remove(self.ParamListSizer)
                self.ParamListSizer = None
            if self.RunButton is not None:
                self.RunButton.Destroy()
                self.RunButton = None
        else:
            self.ParamSizer.Layout()
        self.GetSizer().Layout()
        self.Refresh()
        
    def OnBuildTable(self, event=None):
        names = []
        values = []
        #make names a list of strings
        #make values a list of lists of values
        for param in self.ParamSizer.GetChildren():
            name, val = param.Window.get_values()
            names.append(name)
            values.append(val)
        
        #Build the list of parameters for the parametric study
        if self.ParamListSizer is None:
            #Build and add a sizer for the para values
            self.ParamListBox = wx.StaticBox(self, label = "Parametric Terms Ranges")
            self.ParamListSizer = wx.StaticBoxSizer(self.ParamListBox, wx.VERTICAL)
        else:
            self.ParaList.Destroy()
            self.GetSizer().Remove(self.ParamListSizer)
#            #Build and add a sizer for the para values
            self.ParamListBox = wx.StaticBox(self, label = "Parametric Runs")
            self.ParamListSizer = wx.StaticBoxSizer(self.ParamListBox, wx.VERTICAL)
            
        self.GetSizer().Add(self.ParamListSizer,1,wx.EXPAND)
        self.ParaList = ParametricCheckList(self,names,values)
            
        self.ParamListSizer.Add(self.ParaList,1,wx.EXPAND)
        self.ParaList.SetMinSize((400,-1))
        self.ParamListSizer.Layout()
        self.GetSizer().Layout()
        self.Refresh() 
        
        if self.RunButton is None:
            self.RunButton = wx.Button(self, label='Run Table')
            self.RunButton.Bind(wx.EVT_BUTTON, self.OnRunTable)
            self.ButtonSizer.Add(self.RunButton)
            self.ButtonSizer.Layout()

    def OnRunTable(self, event=None):
        """
        Actually runs the parametric table
        
        This event can only fire if the table is built
        """
        Main = self.GetTopLevelParent()
        sims=[]
        #Column index 1 is the list of parameters
        self.ParaList.GetColumn(1)
        for Irow in range(self.ParaList.GetItemCount()):
            #Loop over all the rows that are checked
            if self.ParaList.IsChecked(Irow):
                
                vals, Names = [], []
                for Icol in range(self.ParaList.GetColumnCount()-1):
                    vals.append(self.ParaList.GetStringItem(Irow, Icol))
                    Names.append(self.ParaList.GetColumn(Icol+1).Text)
                    
                attrs = [self._get_attr(Name) for Name in Names]
                
                # Run the special handler for any additional terms that are
                # not handled in the conventional way using self.items in the 
                # panel.  This is meant for optional terms primarily
                #
                # It can set terms in the GUI so that they can be loaded back by the 
                # simulation builder
                #
                # apply_additional_parametric_terms returns a tuple of attrs, vals 
                # for the terms that were unmatched by the parametric
                # preprocessors
                try:
                    attrs, vals = Main.MTB.InputsTB.apply_additional_parametric_terms(attrs, vals, self.variables)
                except ValueError:
                    break
                
                #Build the recip or the scroll using the GUI parameters
                if Main.SimType == 'recip':
                    sim = Main.build_recip()
                elif Main.SimType == 'scroll':
                    sim = Main.build_scroll()
                else:
                    raise AttributeError
                    
                for val, attr in zip(vals, attrs):
                    # Actually set it
                    setattr(sim, attr, float(val))
                    #Run the post_set_params for all the panels
                    Main.MTB.InputsTB.post_set_params(sim)
                
                #Add an index for the run so that it can be sorted properly
                sim.run_index = Irow + 1
                #Add sim to the list
                sims.append(sim)
        
        #Actually run the batch with the sims that have been built
        Main.run_batch(sims)
        
    def post_prep_for_configfile(self):
        """
        This panel's outputs for the save file
        """
        s = ''
        for i, param in enumerate(self.ParamSizer.GetChildren()):
            name, vals = param.Window.get_values()
            values = ';'.join([str(val) for val in vals])
            s += 'Term' + str(i+1) + ' = Term,' + name +',' + values + '\n'
        return s
    
    def post_get_from_configfile(self,key,value):
        #value is something like "Term1,Piston diameter [m],0.02;0.025"
        string_, value = value.split(',')[1:3]
        #value = Piston diameter [m],0.02;0.025
        #Add a new entry to the table
        self.OnAddTerm()
        I = len(self.ParamSizer.GetChildren())-1
        #Load the values into the variables in the list of variables
        self.ParamSizer.GetItem(I).Window.set_values(string_,value.replace(';',', '))
        
    def _get_attr(self, Name):
        """
        Returns the attribute name corresponding to the given name
        
        Raises
        ------
        ``KeyError`` if not found
        """
        for item in self.variables:
            if item['text'] == Name:
                return item['attr']
        raise KeyError
        
    def update_parametric_terms(self, items):
        """
        Sets the list of possible parametric terms
        """
        self.variables = items
        for child in self.Children:
            if isinstance(child,ParametricOption):
                child.update_parametric_terms(items)

def LabeledItem(parent,id=-1, label='A label', value='0.0', enabled=True, tooltip = None):
    """
    A convenience function that returns a tuple of StaticText and TextCtrl 
    items with the necessary label and values set
    """
    label = wx.StaticText(parent,id,label)
    thing = wx.TextCtrl(parent,id,value)
    if enabled==False:
        thing.Disable()
    if tooltip is not None:
        if enabled:
            thing.SetToolTipString(tooltip)
        else:
            label.SetToolTipString(tooltip)
    return label,thing

class StateChooser(wx.Dialog):
    """
    A dialog used to select the state
    """
    def __init__(self,Fluid,T,rho,parent=None,id=-1,Fluid_fixed = False):
        wx.Dialog.__init__(self,parent,id,"State Chooser",size=(300,250))
        
        class StateChoices(wx.Choicebook):
            def __init__(self, parent, id=-1,):
                wx.Choicebook.__init__(self, parent, id)
                
                self.pageT_dTsh=wx.Panel(self)
                self.AddPage(self.pageT_dTsh,'Saturation Temperature and Superheat')
                self.Tsatlabel1, self.Tsat1 = LabeledItem(self.pageT_dTsh,label="Saturation Temperature [K]",value='290')
                self.Tsat1.default_units = 'Kelvin'
                self.Tsat1.Bind(wx.EVT_CONTEXT_MENU,self.OnChangeUnits)
                self.DTshlabel1, self.DTsh1 = LabeledItem(self.pageT_dTsh,label="Superheat [K]",value='11.1')
                sizer=wx.FlexGridSizer(cols=2,hgap=3,vgap=3)
                sizer.AddMany([self.Tsatlabel1, self.Tsat1,self.DTshlabel1, self.DTsh1])
                self.pageT_dTsh.SetSizer(sizer)
                
                self.pageT_p=wx.Panel(self)
                self.AddPage(self.pageT_p,'Temperature and Absolute Pressure')
                self.Tlabel1, self.T1 = LabeledItem(self.pageT_p,label="Temperature [K]",value='300')
                self.T1.default_units = 'Kelvin'
                self.T1.Bind(wx.EVT_CONTEXT_MENU,self.OnChangeUnits)
                self.plabel1, self.p1 = LabeledItem(self.pageT_p,label="Pressure [kPa]",value='300')
                self.p1.default_units = 'kPa'
                self.p1.Bind(wx.EVT_CONTEXT_MENU,self.OnChangeUnits)
                sizer=wx.FlexGridSizer(cols=2,hgap=3,vgap=3)
                sizer.AddMany([self.Tlabel1, self.T1,self.plabel1, self.p1])
                self.pageT_p.SetSizer(sizer)
            
            def OnChangeUnits(self, event):
                TextCtrl = event.GetEventObject()
                dlg = UnitConvertor(value = float(TextCtrl.GetValue()),
                                    default_units = TextCtrl.default_units
                                    )
                dlg.ShowModal()
                TextCtrl.SetValue(dlg.get_value())
                dlg.Destroy()
        
        sizer=wx.BoxSizer(wx.VERTICAL)
        self.Fluidslabel = wx.StaticText(self,-1,'Fluid: ')
        self.Fluids = wx.ComboBox(self,-1)
        self.Fluids.AppendItems(sorted(CoolProp.__fluids__))
        self.Fluids.SetEditable(False)
        self.Fluids.SetValue(Fluid)
        if Fluid_fixed:
            self.Fluids.Enable(False)
        
        hs = wx.BoxSizer(wx.HORIZONTAL)
        hs.AddMany([self.Fluidslabel,self.Fluids])
        sizer.Add(hs)
        
        sizer.Add((5,5))
        
        self.SC=StateChoices(self)
        sizer.Add(self.SC,1,wx.EXPAND)                
        
        fgs= wx.FlexGridSizer(cols=2,hgap=3,vgap=3)
        self.Tlabel, self.T = LabeledItem(self,label="Temperature [K]",value='300',enabled=False)
        self.plabel, self.p = LabeledItem(self,label="Pressure [kPa]",value='300',enabled=False)
        self.rholabel, self.rho = LabeledItem(self,label="Density [kg/m³]",value='1',enabled=False)
        fgs.AddMany([self.Tlabel,self.T,self.plabel,self.p,self.rholabel,self.rho])
        sizer.Add(fgs)
        
        self.cmdAccept = wx.Button(self,-1,"Accept")
        sizer.Add(self.cmdAccept)
        
        self.SetSizer(sizer)
        self.Fluids.SetStringSelection(Fluid)
        
        if CP.Props(Fluid,"Ttriple") < T < CP.Props(Fluid,"Tcrit"):
            #Pressure from temperature and density
            p = CP.Props('P','T',T,'D',rho,Fluid)
            #Saturation temperature
            Tsat = CP.Props('T','P',p,'Q',1,Fluid)
            self.SC.Tsat1.SetValue(str(Tsat))
            self.SC.DTsh1.SetValue(str(T-Tsat))
            self.SC.T1.SetValue(str(T))
            self.SC.p1.SetValue(str(p))
            self.SC.SetSelection(0) ## The page of Tsat,DTsh
        else:
            #Pressure from temperature and density
            p = CP.Props('P','T',T,'D',rho,Fluid)
            self.SC.T1.SetValue(str(T))
            self.SC.p1.SetValue(str(p))
            self.SC.SetSelection(1) ## The page of Tsat,DTsh
        
        self.OnUpdateVals()
        
        self.SC.Tsat1.Bind(wx.EVT_KEY_UP,self.OnUpdateVals)
        self.SC.DTsh1.Bind(wx.EVT_KEY_UP,self.OnUpdateVals)
        self.SC.T1.Bind(wx.EVT_KEY_UP,self.OnUpdateVals)
        self.SC.p1.Bind(wx.EVT_KEY_UP,self.OnUpdateVals)
        
        self.Fluids.Bind(wx.EVT_COMBOBOX, self.OnFlushVals)
        self.Bind(wx.EVT_CLOSE,self.CancelValues)
        self.cmdAccept.Bind(wx.EVT_BUTTON,self.AcceptValues)
        
        #Bind a key-press event to all objects to get Esc 
        children = self.GetChildren()
        for child in children:
            child.Bind(wx.EVT_KEY_UP,  self.OnKeyPress)
        
    def OnFlushVals(self,event=None):
        """ Clear all the values"""
        self.SC.Tsat1.SetValue("")
        self.SC.DTsh1.SetValue("")
        self.SC.T1.SetValue("")
        self.SC.p1.SetValue("")
        self.T.SetValue("")
        self.p.SetValue("")
        self.rho.SetValue("")
        
    def OnKeyPress(self,event=None):
        """ cancel if Escape key is pressed """
        event.Skip()
        if event.GetKeyCode() == wx.WXK_ESCAPE:
            self.EndModal(wx.ID_CANCEL)
    
    def CancelValues(self,event=None):
        self.EndModal(wx.ID_CANCEL)
        
    def AcceptValues(self,event=None):
        """ If the state is in the vapor phase somewhere, accept it and return """
        Fluid = str(self.Fluids.GetStringSelection())
        T=float(self.T.GetValue())
        p=float(self.p.GetValue())
        if CP.Phase(Fluid,T,p) not in ['Gas','Supercritical']:
            dlg = wx.MessageDialog(None, message = "The phase is not gas or supercritical, cannot accept this state",caption='Invalid state')
            dlg.ShowModal()
            dlg.Destroy()
            return
        self.EndModal(wx.ID_OK)
    
    def GetValues(self):
        Fluid=str(self.Fluids.GetStringSelection())
        T=float(self.T.GetValue())
        p=float(self.p.GetValue())
        rho=float(self.rho.GetValue())
        return Fluid,T,p,rho
        
    def OnUpdateVals(self,event=None):
        if event is not None:
            event.Skip()
            
        PageNum = self.SC.GetSelection()
        Fluid = str(self.Fluids.GetStringSelection())
        try:
            if PageNum == 0:
                #Sat temperature and superheat are given
                p=CP.Props('P','T',float(self.SC.Tsat1.GetValue()),'Q',1.0,Fluid)
                T=float(self.SC.Tsat1.GetValue())+float(self.SC.DTsh1.GetValue())
                rho=CP.Props('D','T',T,'P',p,Fluid)
            elif PageNum == 1:
                #Temperature and pressure are given
                T=float(self.SC.T1.GetValue())
                p=float(self.SC.p1.GetValue())
                rho=CP.Props('D','T',T,'P',p,Fluid)
            else:
                raise NotImplementedError
            
            self.T.SetValue(str(T))
            self.p.SetValue(str(p))
            self.rho.SetValue(str(rho))
        except ValueError:
            return
        

    
class StatePanel(wx.Panel):
    """
    This is a generic Panel that has the ability to select a state given by 
    Fluid, temperature and density by selecting the desired set of inputs in a
    dialog which can be Tsat and DTsh or T & p.
    """
    def __init__(self,parent,id=-1,Fluid='R404A',T=283.15,rho=5.74, Fluid_fixed = False):
        wx.Panel.__init__(self,parent,id)
        
        self._Fluid_fixed = Fluid_fixed
        p = CP.Props('P','T',T,'D',rho,str(Fluid))
        sizer=wx.FlexGridSizer(cols=2,hgap=4,vgap=4)
        self.Fluidlabel, self.Fluid = LabeledItem(self,label="Fluid",value=str(Fluid))
        
        self.Tlabel, self.T = LabeledItem(self,label="Temperature [K]",value=str(T))
        self.plabel, self.p = LabeledItem(self,label="Pressure [kPa]",value=str(p))
        self.rholabel, self.rho = LabeledItem(self,label="Density [kg/m³]",value=str(rho))
        sizer.AddMany([self.Fluidlabel, self.Fluid,self.Tlabel,self.T,self.plabel,self.p,self.rholabel,self.rho])
        
        
        for box in [self.T,self.p,self.rho,self.Fluid]:
            #Make the box not editable
            self.Fluid.SetEditable(False)
            #Bind events tp fire the chooser when text boxes are clicked on
            box.Bind(wx.EVT_LEFT_DOWN,self.UseChooser)
            box.SetToolTipString('Click on me to select the state')
        
        self.SetSizer(sizer)
        
    def GetState(self):
        """
        returns a :class:`CoolProp.State.State` instance from the given values
        in the panel
        """
        Fluid = str(self.Fluid.GetValue())
        T = float(self.T.GetValue())
        rho = float(self.rho.GetValue())
        return State(Fluid,dict(T=T,D=rho))
    
    def set_state(self, Fluid, **kwargs):
        """
        Fluid must not be unicode
        """
        if self._Fluid_fixed and not str(self.Fluid.GetValue()) == str(Fluid):
            import warnings
            warnings.warn('Could not set state since fluid is fixed')
            return

        #Create a state instance from the parameters passed in
        S  = State(str(Fluid), kwargs)

        #Load up the textboxes
        self.Fluid.SetValue(S.Fluid)
        self.T.SetValue(str(S.T))
        self.rho.SetValue(str(S.rho))
        self.p.SetValue(str(S.p))
    
    def UseChooser(self,event=None):
        """
        An event handler that runs the State Chooser dialog and sets the
        values back in the panel
        """
        #Values from the GUI
        Fluid = str(self.Fluid.GetValue())
        T = float(self.T.GetValue())
        rho = float(self.rho.GetValue())
        
        #Instantiate the chooser Dialog
        SCfrm=StateChooser(Fluid=Fluid,T=T,rho=rho,Fluid_fixed = self._Fluid_fixed)
        
        #If they clicked accept
        if wx.ID_OK == SCfrm.ShowModal():
            Fluid,T,p,rho=SCfrm.GetValues()
            #Set the GUI values
            self.Fluid.SetValue(str(Fluid))
            self.T.SetValue(str(T))
            self.p.SetValue(str(p))
            self.rho.SetValue(str(rho))
        SCfrm.Destroy()

class StateInputsPanel(PDPanel):
    
    def __init__(self, parent, configfile, **kwargs):
    
        PDPanel.__init__(self, parent, **kwargs)
        
        #Loads all the parameters from the config file (case-sensitive)
        self.configdict, self.descdict = self.get_from_configfile('StatePanel')
        
        self.items = [
                      dict(attr='omega')
                      ]
        
        box_sizer = wx.BoxSizer(wx.VERTICAL)
        
        sizer = wx.FlexGridSizer(cols=2, vgap=4, hgap=4)
        self.ConstructItems([self.items[0]],sizer,self.configdict,self.descdict)
        box_sizer.Add(sizer)
        
        box_sizer.Add(wx.StaticText(self,-1,"Suction State"))
        box_sizer.Add(wx.StaticLine(self, -1, (25, 50), (300,1)))    
            
        Fluid = self.configdict['inletState']['Fluid']
        T = self.configdict['inletState']['T']
        rho = self.configdict['inletState']['rho']
        self.SuctionState = StatePanel(self,Fluid=Fluid,T=T,rho=rho)
        box_sizer.Add(self.SuctionState)
        
        box_sizer.Add((20,20))
        box_sizer.Add(wx.StaticText(self,-1,"Discharge State"))
        box_sizer.Add(wx.StaticLine(self,-1,(25, 50), (300,1)))
        
        self.cmbDischarge = wx.ComboBox(self)
        self.cmbDischarge.AppendItems(['Discharge pressure [kPa]', 
                                       'Pressure ratio [-]', 
                                       'Saturated Temperature [K]'
                                       ])
        self.cmbDischarge.SetStringSelection(self.Discharge_key)
        self.DischargeValue = wx.TextCtrl(self, value = self.Discharge_value)
        self.cmbDischarge.Bind(wx.EVT_COMBOBOX, self.OnChangeDischargeVariable)
        self.DischargeValue.Bind(wx.EVT_KILL_FOCUS,self.OnChangeDischargeValue)
        self.DischargeValue.Bind(wx.EVT_TEXT,self.OnChangeDischargeValue)
        
        sizer = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
        sizer.AddMany([self.cmbDischarge, self.DischargeValue])
        box_sizer.Add(sizer)
        
        self.SetSizer(box_sizer)
        sizer.Layout()  
        
        #Set the value self._discharge_pressure variable
        self.OnChangeDischargeValue()
        
    def OnChangeDischargeValue(self, event = None):
        """ 
        Set the internal pressure variable when the value is changed in the TextCtrl
        """
        p_suction = self.SuctionState.GetState().p
        Fluid = self.SuctionState.GetState().Fluid
        
        if self.cmbDischarge.GetStringSelection() == 'Discharge pressure [kPa]':
            p_disc = float(self.DischargeValue.GetValue())
            self._discharge_pressure = p_disc
            
        elif self.cmbDischarge.GetStringSelection() == 'Pressure ratio [-]':
            p_ratio = float(self.DischargeValue.GetValue())
            self._discharge_pressure = p_ratio*p_suction
        
        elif self.cmbDischarge.GetStringSelection() == 'Saturated Temperature [K]':
            Tsat = float(self.DischargeValue.GetValue())
            self._discharge_pressure = CP.Props('P','T',Tsat,'Q',1.0,Fluid)

        else:
            raise KeyError
        
    def OnChangeDischargeVariable(self, event = None):
        """
        Fires when the combobox is selected
        """
        p_suction = self.SuctionState.GetState().p
        Fluid = self.SuctionState.GetState().Fluid

        #Remove the handler for units selection
        self.DischargeValue.Unbind(wx.EVT_CONTEXT_MENU)
        
        if self.cmbDischarge.GetStringSelection() == 'Discharge pressure [kPa]':
            self.DischargeValue.SetValue(str(self._discharge_pressure))
            #Set the handler for unit selection
            self.DischargeValue.default_units = 'kPa'
            self.DischargeValue.Bind(wx.EVT_CONTEXT_MENU, self.OnChangeUnits)
            
        elif self.cmbDischarge.GetStringSelection() == 'Pressure ratio [-]':
            p_disc = self._discharge_pressure
            pratio = p_disc/p_suction
            self.DischargeValue.SetValue(str(pratio))
        
        elif self.cmbDischarge.GetStringSelection() == 'Saturated Temperature [K]':
            p_disc = self._discharge_pressure
            Tsat = CP.Props('T', 'P', p_disc, 'Q', 1.0, Fluid)
            self.DischargeValue.SetValue(str(Tsat))
            #Set the handler for unit selection
            self.DischargeValue.default_units = 'Kelvin'
            self.DischargeValue.Bind(wx.EVT_CONTEXT_MENU, self.OnChangeUnits)
        
        else:
            raise KeyError
        
    def post_get_from_configfile(self, key, value):
        Dummy, value, key = value.split(',')
        self.Discharge_key = key
        self.Discharge_value = str(value)
        
    def post_set_params(self, simulation):
        Fluid = self.SuctionState.GetState().Fluid
        
        simulation.inletState = self.SuctionState.GetState()
            
        simulation.discharge_pressure = self._discharge_pressure
            
        #Set the state variables in the simulation
        simulation.suction_pressure = self.SuctionState.GetState().p
        simulation.suction_temp = self.SuctionState.GetState().T
        
        if self.SuctionState.GetState().p < CP.Props(Fluid,'pcrit'):
            #If subcritical, also calculate the superheat and sat_temp
            p = simulation.suction_pressure
            simulation.suction_sat_temp = CP.Props('T', 'P', p, 'Q', 1.0, Fluid)
            simulation.suction_superheat = simulation.suction_temp-simulation.suction_sat_temp
        else:
            #Otherwise remove the parameters
            del simulation.suction_sat_temp
            del simulation.suction_superheat
            
        #Set the state variables in the simulation
        simulation.discharge_pratio = simulation.discharge_pressure/simulation.suction_pressure
        
        if simulation.discharge_pressure < CP.Props(Fluid,'pcrit'):
            p = simulation.discharge_pressure
            simulation.discharge_sat_temp = CP.Props('T', 'P', p, 'Q', 1.0, Fluid)
        else:
            del simulation.discharge_sat_temp
        
    def post_prep_for_configfile(self):
        """
        Write a string representation of the state
        """
        State_ = self.SuctionState.GetState()
        StateString = 'inletState = State,'+State_.Fluid+','+str(State_.T)+','+str(State_.rho)
        DischargeString = 'discharge = Discharge,'+str(self.DischargeValue.GetValue())+','+self.cmbDischarge.GetStringSelection()
        return StateString+'\n'+DischargeString+'\n'
    
    def get_additional_parametric_terms(self):
        return [dict(attr = 'suction_pressure',
                     text = 'Suction pressure [kPa]',
                     parent = self),
                dict(attr = 'suction_sat_temp',
                     text = 'Suction saturated temperature (dew) [K]',
                     parent = self),
                dict(attr = 'suction_temp',
                     text = 'Suction temperature [K]',
                     parent = self),
                dict(attr = 'suction_superheat',
                     text = 'Superheat [K]',
                     parent = self),
                dict(attr = 'discharge_pressure',
                     text = 'Discharge pressure [kPa]',
                     parent = self),
                dict(attr = 'discharge_sat_temp',
                     text = 'Discharge saturated temperature (dew) [K]',
                     parent = self),
                dict(attr = 'discharge_pratio',
                     text = 'Discharge pressure ratio [-]',
                     parent = self)
                ]
    
    def apply_additional_parametric_terms(self, attrs, vals, panel_items):
        """
        Set parametric terms for the state panel based on parameters obtained
        from the parametric table
        """
        panel_attrs = [panel_item['attr'] for panel_item in panel_items]
        # First check about the suction state; if two suction related terms are 
        # provided, use them to fix the inlet state
        suct_params = [(par,val) for par,val in zip(attrs,vals) if par.startswith('suction')]
        num_suct_params = len(suct_params)
        
        #Get a copy of the state from the StatePanel
        inletState = self.SuctionState.GetState()
        
        if num_suct_params>0:
            #Unzip the parameters (List of tuples -> tuple of lists)
            suct_attrs, suct_vals = zip(*suct_params)
            
        if num_suct_params == 2:
            # Remove all the entries that correspond to the suction state - 
            # we need them and don't want to set them in the conventional way
            for a in suct_attrs:
                i = attrs.index(a)
                vals.pop(i)
                attrs.pop(i)
            
            #Temperature and pressure provided
            if 'suction_temp' in suct_attrs and 'suction_pressure' in suct_attrs:
                suction_temp = suct_vals[suct_attrs.index('suction_temp')]
                suction_pressure = suct_vals[suct_attrs.index('suction_pressure')]
                self.SuctionState.set_state(inletState.Fluid,
                                            T=suction_temp, 
                                            P=suction_pressure)
                
            #Dew temperature and superheat provided
            elif 'suction_sat_temp' in suct_attrs and 'suction_superheat' in suct_attrs:
                suction_sat_temp = suct_vals[suct_attrs.index('suction_sat_temp')]
                suction_superheat = suct_vals[suct_attrs.index('suction_superheat')]
                suction_temp = suction_sat_temp + suction_superheat
                suction_pressure = CP.Props('P','T',suction_sat_temp,'Q',1.0,inletState.Fluid)
                self.SuctionState.set_state(inletState.Fluid,
                                            T=suction_temp, 
                                            P=suction_pressure)
            else:
                raise ValueError('Invalid combination of suction states: '+str(suct_attrs))
            
        elif num_suct_params == 1:
            raise NotImplementedError('only one param provided')
        elif num_suct_params >2:
            raise ValueError ('Only two inlet state parameters can be provided in parametric table')
        
        # Then check about the discharge state; only one variable is allowed
        disc_params = [(par,val) for par,val in zip(attrs,vals) if par.startswith('discharge')]
        num_disc_params = len(disc_params)
        
        if num_disc_params>0:
            #Unzip the parameters (List of tuples -> tuple of lists)
            disc_attrs, disc_vals = zip(*disc_params)
            
        if num_disc_params == 1:
            # Remove all the entries that correspond to the discharge state - 
            # we need them and don't want to set them in the conventional way
            for a in disc_attrs:
                i = attrs.index(a)
                vals.pop(i)
                attrs.pop(i)
                
            if 'discharge_pressure' in disc_attrs:
                discharge_pressure = disc_vals[disc_attrs.index('discharge_pressure')]
                self._discharge_pressure = discharge_pressure
                
            elif 'discharge_pratio' in disc_attrs:
                discharge_pratio = disc_vals[disc_attrs.index('discharge_pratio')]
                p_suction = self.SuctionState.GetState().p
                self._discharge_pressure = discharge_pratio * p_suction
        
            elif 'discharge_sat_temp' in disc_attrs:
                Tsat = disc_vals[disc_attrs.index('discharge_sat_temp')]
                Fluid = inletState.Fluid
                self._discharge_pressure = CP.Props('P','T',Tsat,'Q',1.0,Fluid)
            
            #Fire the event manually to update the textbox
            self.OnChangeDischargeVariable()
            
        elif num_disc_params > 1:
            raise ValueError ('Only one discharge pressure parameter can be provided in parametric table')
            
        return attrs, vals
        
if __name__=='__main__':
    execfile('PDSimGUI.py')