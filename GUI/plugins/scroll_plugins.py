import sys
sys.path.append('..')

import pdsim_plugins
from panels import pdsim_panels
from PDSim.flow.flow import FlowPath
from PDSim.core.containers import Tube
from PDSim.core.containers import ControlVolume
from PDSim.scroll import scroll_geo
from math import pi
import wx
from wx.lib.scrolledpanel import ScrolledPanel
from PDSim.scroll.plots import plotScrollSet, ScrollAnimForm
import types
import warnings

LabeledItem = pdsim_panels.LabeledItem

class struct(object):
    pass
                
class InjectionPortPanel(wx.Panel):
    """
    A panel with the values for one injection port
    """
    def __init__(self, parent, index):
        wx.Panel.__init__(self,parent)
        
        #: The index of the panel
        self.index = index
        
        #: A textual string with the 
        self.indexText = wx.StaticText(self,label='#'+str(index))
        
        self.AddPort = wx.Button(self, label='+',style = wx.ID_REMOVE)
        self.RemovePort = wx.Button(self,label='-',style = wx.ID_REMOVE)
        self.AddPort.Bind(wx.EVT_BUTTON,lambda(event):self.Parent.OnAddPort(self))
        self.RemovePort.Bind(wx.EVT_BUTTON,lambda(event):self.Parent.OnRemovePort(self))
        
        button_sizer = wx.BoxSizer(wx.VERTICAL)
        button_sizer.Add(self.indexText,0,wx.EXPAND| wx.ALIGN_CENTER)
        button_sizer.Add(self.AddPort,0,wx.EXPAND| wx.ALIGN_CENTER)
        button_sizer.Add(self.RemovePort,0,wx.EXPAND| wx.ALIGN_CENTER)
        
        element_sizer = wx.FlexGridSizer(cols=2,hgap = 2, vgap = 2)
        element_sizer.AddSpacer(10)
        element_sizer.AddSpacer(10)
        element_sizer.Add(wx.StaticText(self,label="Diameter [m]"))
        self.D_port = wx.TextCtrl(self,value="0.001") 
        self.D_port.SetToolTipString('Diameter of the port')
        element_sizer.Add(self.D_port)
        element_sizer.Add(wx.StaticText(self,label="Offset [m]"))
        self.offset_port = wx.TextCtrl(self,value="0.001") 
        self.offset_port.SetToolTipString('Offset of the port away from involute it is referenced from')
        element_sizer.Add(self.offset_port)
        element_sizer.Add(wx.StaticText(self,label="Involute Angle"))
        self.phi_inj_port = wx.TextCtrl(self,value="7.141") 
        self.phi_inj_port.SetToolTipString('If you want symmetric injection ports, the involute angle on the inner involute should be pi radians greater than that on the outer involute.  View the ports to be sure')
        element_sizer.Add(self.phi_inj_port)
        element_sizer.Add(wx.StaticText(self,label="Neighbor Involute"))
        self.involute = wx.Choice(self)
        self.involute.AppendItems(['Outer involute','Inner involute'])
        self.involute.SetSelection(0)
        element_sizer.Add(self.involute)
        element_sizer.Add(wx.StaticText(self,label="X_d"))
        self.X_d = wx.TextCtrl(self,value="0.8") 
        self.X_d.SetToolTipString('Flow coefficient into the working chamber')
        element_sizer.Add(self.X_d)
        element_sizer.Add(wx.StaticText(self,label="X_d_backflow"))
        self.X_d_backflow = wx.TextCtrl(self,value="0.8") 
        self.X_d_backflow.SetToolTipString('Flow coefficient for backflow')
        element_sizer.Add(self.X_d_backflow)
        
        element_sizer.Add(wx.StaticText(self,label="Uses check valve"))
        self.check_valve = wx.CheckBox(self,label="")
        element_sizer.Add(self.check_valve)
        self.SymmLabel = wx.StaticText(self,label='Symmetric with #')
        self.SymmTarget = wx.ComboBox(self)
        self.SymmTarget.Append('None')
        self.SymmTarget.SetSelection(0)
        self.SymmTarget.SetEditable(False)
        element_sizer.Add(self.SymmLabel)
        element_sizer.Add(self.SymmTarget)
        self.SymmTarget.Bind(wx.EVT_COMBOBOX, self.OnMakeSymmetric)
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(button_sizer,1,wx.ALIGN_CENTER_VERTICAL)
        sizer.Add(element_sizer)
        self.SetSizer(sizer)
        sizer.Layout()
        
        self.set_index(index)
        
    def set_index(self, index):
        """
        A convenience function for setting the index for the port and 
        doing any further work required in the GUI
        
        Parameters
        ----------
        index : int
            The 1-based index of the port
        """
        self.index = index
        self.indexText.SetLabel('#'+str(index))
        
    def set_values(self,vals):
        """
        Takes in a dictionary of terms and sets the values in the panel
        
        Parameters
        ----------
        phi : float
        involute : string
            If ``'i'``, phi is along the inner involute of the fixed scroll
            If ``'o'``, phi is along the outer involute of the fixed scroll
        check_valve : boolean
        symmetric : string
            If 'None', no symmetric link for this port, otherwise it is the string
            representation of an integer with the 1-based index of the port
        D : float
            Diameter of the port
        offset : float
            Offset from the scroll wrap
        """
        self.phi_inj_port.SetValue(str(vals['phi']))
        self.check_valve.SetValue(vals['check_valve'])
        self.D_port.SetValue(str(vals['D']))
        self.offset_port.SetValue(str(vals['offset']))
        if vals['involute'] == 'i':
            self.involute.SetStringSelection('Inner involute')
        elif vals['involute'] == 'o':
            self.involute.SetStringSelection('Outer involute')
        else:
            raise ValueError
        self.SymmTarget.SetStringSelection(str(vals['symmetric']))
        if not vals['symmetric'] == 'None':
            self.OnMakeSymmetric()
        self.X_d.SetValue(str(vals['X_d']))
        self.X_d_backflow.SetValue(str(vals['X_d_backflow']))
    
    def get_values(self):
        """
        Returns a dictionary of values
        
        Variables as described as in set_values()
        """
        if self.involute.GetStringSelection() == 'Outer involute':
            inner_outer = 'o'
        elif self.involute.GetStringSelection() == 'Inner involute':
            inner_outer = 'i'
        else:
            raise ValueError
        
        return dict(phi = float(self.phi_inj_port.GetValue()),
                    involute = inner_outer,
                    symmetric = self.SymmTarget.GetStringSelection(),
                    check_valve = self.check_valve.IsChecked(),
                    D = float(self.D_port.GetValue()),
                    offset = float(self.offset_port.GetValue()),
                    X_d = float(self.X_d.GetValue()),
                    X_d_backflow = float(self.X_d_backflow.GetValue())
                    )

    
    def OnMakeSymmetric(self, event = None):
        """
        An event handler for changing the symmetric nature of the 
        """
        # Skip the event so that other objects also running the same handler can
        # do so for instance if there are multiple ports that are symmetric with
        # the port
        if event is not None:
            event.Skip()
        #Get the index if possible 
        I = None   
        try:
            I = int(self.SymmTarget.GetStringSelection())
        except ValueError:
            pass
        
        #If symmetric linking is enabled, 
        if I is not None:
            #Get the symmetric port and its involute angle
            port = self.Parent.ports_list[I-1]
            phi = float(port.phi_inj_port.GetValue())
            if port.involute.GetStringSelection() == 'Inner involute':
                self.involute.SetStringSelection('Outer involute')
                self.phi_inj_port.SetValue(str(phi-pi))
            if port.involute.GetStringSelection() == 'Outer involute':
                self.involute.SetStringSelection('Inner involute')
                self.phi_inj_port.SetValue(str(phi+pi))
            self.involute.Enable(False)
            self.phi_inj_port.Enable(False)
            port.SymmTarget.Enable(False)
            
            #If this port instance fired the event (ie not from the symmetric port)
            #bind the events to the symmetric object
            if event is None or event.GetEventObject().Parent == self:
                port.involute.Bind(wx.EVT_COMBOBOX, self.OnMakeSymmetric)
                port.phi_inj_port.Bind(wx.EVT_TEXT, self.OnMakeSymmetric)
        else:
            self.involute.Enable(True)
            self.phi_inj_port.Enable(True)
        
class InjectionElementPanel(wx.Panel):
    """
    A panel with the injection values for one injection line with a box around it
    """
    def __init__(self, parent,index):
        wx.Panel.__init__(self,parent)
        
        #Inputs Toolbook
        ITB = self.GetTopLevelParent().MTB.InputsTB
        CPState = None
        for panel in ITB.panels:
            if panel.Name == 'StatePanel':
                CPState = panel.SuctionStatePanel.GetState()
                break
        if CPState is None:
            raise ValueError('StatePanel not found in Inputs Toolbook')
        
        #You can only inject the same refrigerant as at the suction so fix the fluid
        self.state = pdsim_panels.StatePanel(self, CPState=CPState, Fluid_fixed = True)
        
        self.Llabel,self.Lval = LabeledItem(self, label='Length of injection line',value='1.0')
        self.IDlabel,self.IDval = LabeledItem(self, label='Inner diameter of injection line',value='0.01')
        line_sizer = wx.FlexGridSizer(cols = 2)
        line_sizer.AddMany([self.Llabel,self.Lval])
        line_sizer.AddMany([self.IDlabel,self.IDval])
        
        self.RemoveButton = wx.Button(self,label='Remove Line')
        #Parent of IEP is scrolled_panel, need to go up to IIP
        self.RemoveButton.Bind(wx.EVT_BUTTON, lambda event: self.GrandParent.RemoveInjection(self))
        
        self.SizerBox = wx.StaticBox(self, label = "Injection line #"+str(index))
        self.SBSSizer = wx.StaticBoxSizer(self.SizerBox, wx.VERTICAL)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.RemoveButton)
        
        text = wx.StaticText(self, label='Parameters:')
        font = text.GetFont()
        font.SetUnderlined(True)
        text.SetFont(font)
        sizer.Add(text)
        sizer.AddSpacer(10)
        sizer.Add(line_sizer)
        sizer.AddSpacer(10)
        text = wx.StaticText(self, label='State:')
        font = text.GetFont()
        font.SetUnderlined(True)
        text.SetFont(font)
        sizer.Add(text)
        sizer.Add(self.state)
        text = wx.StaticText(self, label='Injection Ports:')
        font = text.GetFont()
        font.SetUnderlined(True)
        text.SetFont(font)
        sizer.Add(text)
        IPP = InjectionPortPanel(self, index = 1)
        self.ports_list = [IPP]
        sizer.Add(IPP)
        self.InnerSizer = sizer
        self.SBSSizer.Add(sizer)
        self.SBSSizer.Layout()
        self.SetSizer(self.SBSSizer)
        self.Nports = 1
        
    def OnRemovePort(self, removed_port = None):
        #Store a dictionary of partners for each port
        partner_dict = {}
        for port in self.ports_list:
            #Get the index
            index_partner = port.SymmTarget.GetStringSelection()
            try:
                #Get the partner panel
                partner = self.ports_list[int(index_partner)-1]
            except ValueError:
                partner = None
                
            #If partner still exists after removing the port, keep a reference
            #to the port
            if partner == removed_port:
                #Entry points to nothing
                partner_dict[port.index] = None
            else:
                #Entry points to partner instance
                partner_dict[port.index] = partner
        
        #Can't remove the first element
        if self.Nports > 1:
            I_removed = self.ports_list.index(removed_port)
            self.ports_list.remove(removed_port)
            removed_port.Destroy()
            self.SBSSizer.Layout()
            self.Refresh()
            self.Nports -= 1
            self.Fit()
            
            #renumber the ports starting at 1
            for i,port in enumerate(self.ports_list):
                port.set_index(i+1)
            
            # Shift the indices down for all the elements above the removed port
            # in the partner_dictionary
            partner_dict.pop(I_removed+1) #partner_dict uses 1-based indexing
            for i in range(I_removed+1,len(self.ports_list)+1):
                partner_dict[i] = partner_dict[i+1]
            partner_dict.pop(len(self.ports_list)+1)
                
            Nports = len(self.ports_list)
            for port in self.ports_list:
                    
                port.SymmTarget.Clear()
                port.SymmTarget.Append('None')
                for i in range(1, Nports+1):
                    if not port.index == i:
                        #Add this new port to the list of possible partners
                        port.SymmTarget.Append(str(i))
                
                partner = partner_dict[port.index]
                #Reset the partner if it had a partner
                if partner is not None:
                    port.SymmTarget.SetStringSelection(str(partner.index))
                else:
                    port.SymmTarget.SetStringSelection('None')
                    
            return True
        else:
            return False
        
        self.Parent.FitInside()
        
    def OnAddPort(self, event = None):
        IPP = InjectionPortPanel(self, index = len(self.ports_list)+1)
        self.SBSSizer.Add(IPP)
        self.ports_list.append(IPP)
        self.SBSSizer.Layout()
        self.Nports += 1
        self.GetSizer().Layout()
        self.Fit()
        self.Refresh()
        
        # When you add a port, it cannot be a partner of any other chamber at
        # instantiation
        Nports = len(self.ports_list)
        for port in self.ports_list:
            
            index_partner = None
            if not port.SymmTarget.GetStringSelection() == 'None':
                #Get the index
                index_partner = port.SymmTarget.GetStringSelection()
                
            port.SymmTarget.Clear()
            port.SymmTarget.Append('None')
            for i in range(1, Nports+1):
                if not port.index == i:
                    #Add this new port to the list of possible partners
                    port.SymmTarget.Append(str(i))
            
            #Reset the value if it had a partner
            if index_partner is not None:
                port.SymmTarget.SetStringSelection(index_partner)
            else:
                port.SymmTarget.SetStringSelection('None')
        
        self.Parent.FitInside()
        
class InjectionInputsPanel(pdsim_panels.PDPanel):
    """
    The container panel for all the injection ports and injection data 
    """ 
    def __init__(self, parent, **kwargs):
        pdsim_panels.PDPanel.__init__(self,parent,**kwargs)
        
        #Now we are going to put everything into a scrolled window
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        
        self.scrolled_panel = ScrolledPanel(self, size=(-1,-1),
                                 style = wx.TAB_TRAVERSAL, name="panel1")
        self.scrolled_panel.SetScrollbars(1,1,1,1)
        self.scrolled_panel.SetupScrolling()
        
        #Add the header row of buttons
        self.View = wx.Button(self.scrolled_panel, label='View')
        self.View.Bind(wx.EVT_BUTTON, self.OnView)
        self.AddInjection = wx.Button(self.scrolled_panel, label='Add Injection Line')
        self.AddInjection.Bind(wx.EVT_BUTTON, self.OnAddInjection)
        self.PlotExistence = wx.Button(self.scrolled_panel, label='Plot Existence')
        self.PlotExistence.Bind(wx.EVT_BUTTON, self.OnPlotExistence)
        buttons_sizer = wx.BoxSizer(wx.HORIZONTAL)
        buttons_sizer.Add(self.AddInjection)
        buttons_sizer.Add(self.View)
        buttons_sizer.Add(self.PlotExistence)
        
        sizer = wx.FlexGridSizer(cols = 1)
        sizer.Add(buttons_sizer)
        sizer.AddSpacer(10)
        sizer.Layout()
        
        self.scrolled_panel.SetAutoLayout(1)

        #Do the layout of all the panels
        self.scrolled_panel.SetSizer(sizer)
        main_sizer.Add(self.scrolled_panel,1,wx.EXPAND)
        self.SetSizer(main_sizer)
        main_sizer.Layout()
        
        #Set some local variables
        self.Nterms = 0
        self.Lines = []
        
    def OnAddInjection(self, event = None):
        """
        Add an injection line to the injection panel
        """
        IE = InjectionElementPanel(self.scrolled_panel,self.Nterms+1)
        #Put the panel within the scrolled panel and refresh
        self.scrolled_panel.GetSizer().Add(IE,0)
        self.scrolled_panel.FitInside()
        self.GetSizer().Layout()
        
        #Update the local variables
        self.Lines.append(IE)
        self.Nterms += 1
        
        self.Refresh()
    
    def remove_all(self):
        while self.Lines:
            self.RemoveInjection(self.Lines[0])
        
    def RemoveInjection(self, injection):
        """
        Remove the given injection term
        """
        self.Lines.remove(injection)
        injection.Destroy()
        self.Nterms -= 1
        #Renumber the injection panels that are contained in scrolled_panel
        I=1
        for child in self.scrolled_panel.Children:
            if isinstance(child,InjectionElementPanel):
                child.SizerBox.SetLabel("Injection line #"+str(I))
                I+=1
        self.GetSizer().Layout()
        self.scrolled_panel.FitInside()
        self.Refresh()
        
    def OnView(self, event):
        
        geo = self.GetTopLevelParent().MTB.InputsTB.panels[0].Scroll.geo
        SAF = ScrollAnimForm(geo, start = False)
        
        #IEPs are children that are instances of InjectionElementPanel class
        IEPs = [child for child in self.scrolled_panel.Children if isinstance(child, InjectionElementPanel)]
        for IEP in IEPs:
            for child in IEP.Children:
                if isinstance(child, InjectionPortPanel):
                    # Get the values from the panel
                    vals = child.get_values()
                    # Overlay the port on the scroll wrap plot
                    scroll_geo.overlay_injection_port(0, geo, vals['phi'], SAF.ax, vals['involute'], rport = vals['D']/2, offset = vals['offset'])
        SAF.start()
        SAF.Show()
        
    def OnPlotExistence(self, event = None):
        """
        Plot a 2D line plot showing which control volume is connected to 
        each injection port as a function of the crank angle
        """
        import pylab
        import numpy as np
        
        _Scroll = self.GetTopLevelParent().MTB.InputsTB.panels[0].Scroll
        
        Iport = 1
        #IEPs are children that are instances of InjectionElementPanel class
        IEPs = [child for child in self.scrolled_panel.Children if isinstance(child,InjectionElementPanel)]
        for IEP in IEPs:
            for child in IEP.Children:
                if isinstance(child,InjectionPortPanel):
                    #Get the values from the port panel
                    v = child.get_values()
                    
                    partner_list = []
                    
                    theta = np.linspace(0, 2*pi, 1000)
                    for th in theta:
                        partner_list.append(_Scroll._get_injection_CVkey(v['phi'], th, v['involute']))
        
                    #Find the break points in each segment
                    dividers = [i for i in range(len(theta)-1) if not partner_list[i] == partner_list[i+1]]
                    #Add end and beginning indices
                    dividers = [0]+dividers+[len(theta)-1]
                        
                    for i in range(len(dividers)-1):
                        L = dividers[i]
                        R = dividers[i+1]
                        M = int((L + R)/2)
                        pylab.plot(np.r_[theta[L],theta[R]],np.r_[Iport,Iport])
                        pylab.plot(np.r_[theta[L],theta[L]],np.r_[Iport-0.02,Iport+0.02],'k')
                        pylab.plot(np.r_[theta[R],theta[R]],np.r_[Iport-0.02,Iport+0.02],'k')
                        pylab.text(theta[M],Iport+.02,partner_list[M],ha = 'center', va='bottom')    
                    
                    #Increase the counter
                    Iport += 1
        
        pylab.xticks([0,pi/2,pi,3*pi/2,2*pi],
                     [0,r'$\pi/2$',r'$\pi$',r'$3\pi/2$',r'$2\pi$'])
        pylab.xlim(0,2*pi)
        pylab.ylim(0.5,Iport-1+0.5)
        pylab.yticks(range(1,Iport+1))
        pylab.show()
        
    def build_from_configfile(self, config):
        """
        Get parameters from the configfile section for this plugin
        
        Parameters
        ----------
        config : yaml configuration section for the plugin
        
        """

        if config:
            self.remove_all()
            for line in config:
                # Add an injection line panel
                self.OnAddInjection()
                #Get a pointer to the last IEP (the one just added)
                IEP = self.Lines[-1]
                #Set the line length in the GUI [m]
                IEP.Lval.SetValue(str(line['Length']))
                #Set the line ID in the GUI
                IEP.IDval.SetValue(str(line['ID']))
                #Set the State in the GUI
                State = line['inletState']
                IEP.state.set_state(State['Fluid'], T = State['T'], D = State['rho'])
                if 'ports' in line and line['ports']:
                    for i,port in enumerate(line['ports']):
                        if i > 0: IEP.OnAddPort()
                        # Get a pointer to the port panel
                        portpanel = IEP.ports_list[-1]
                        # Set the values in the panel
                        portpanel.set_values(port)
    
    def get_additional_parametric_terms(self):
        
        #: the list of terms
        _T = []
        
        #IEPs are children of injection_panel that are instances of InjectionElementPanel class
        IEPs = [child for child in self.scrolled_panel.Children if isinstance(child,InjectionElementPanel)]
        for i,IEP in enumerate(IEPs):
            I = str(i+1)
            
            _T += [dict(attr = 'injection_state_pressure_' + I,
                        text = 'Injection pressure #' + I + ' [kPa]',
                        parent = self),
                   dict(attr = 'injection_state_sat_temp_' + I,
                        text = 'Injection saturated temperature (dew) #' + I + ' [K]',
                        parent = self),
                   dict(attr = 'injection_state_temp_' + I,
                        text = 'Injection temperature #' + I + ' [K]',
                        parent = self),
                   dict(attr = 'injection_state_superheat_' + I,
                        text = 'Injection superheat #' + I + ' [K]',
                        parent = self),
                ]
                
            Ports = [c for c in IEP.Children if isinstance(c,InjectionPortPanel)]
            for j,child in enumerate(Ports):
                J = str(j+1)
                _T += [dict(attr = 'injection_phi_'+I+'_'+J,
                            text = 'Injection port angle #'+I+':'+J+' [rad]',
                            parent = self)]
                
        return _T
                
        
    def apply_additional_parametric_terms(self, attrs, vals, panel_items):
        """
        Set the terms in the injection panel based on the additional parametric
        terms provided by the get_additional_parametric_terms() function
        """
        
        def apply_line_terms(attrs, vals):
            
            def is_int(i):
                """ Returns True if it is an integer """
                try:
                    i = int(i)
                    return True
                except ValueError:
                    return False
                
            def is_line_term(attr):
                """
                Check if it is a line type term of the form injection_xxxxx_1'
                and is not a port term of the form injection_xxxxx_1_1
                """
                if not attr.startswith('injection'):
                    return False
                
                #If there are no underscores, return false
                if len(attr.rsplit('_',1)) == 1:
                    return False
                
                #Try to split twice
                attr,i,j = attr.rsplit('_',2)
                
                # If the far right one is an integer and the left part isn't you are
                # ok, its an injection line
                if not is_int(i) and is_int(j):
                    return True
                else:
                    return False
        
            # First check about the injection state; if two state related terms are 
            # provided, use them to fix the injection state
            inj_state_params = [(par,val) for par,val in zip(attrs,vals) if is_line_term(par)]
            num_inj_state_params = len(inj_state_params)
            
            for i in range(len(self.Lines)):
                
                #Find the injection state terms that apply for this line
                state_params = [(par,val) for par,val in zip(attrs,vals) 
                                if par.find('state') > -1 and par.endswith(str(i+1))]
                num_state_params = len(state_params)
                
                #Get a copy of the state from the StatePanel
                inletState = self.Lines[i].state.GetState()
                
                if num_state_params > 0:
                    #Unzip the parameters (List of tuples -> tuple of lists)
                    state_attrs, state_vals = zip(*state_params)
                    
                if num_state_params == 2:
                    # Remove all the entries that correspond to the injection state - 
                    # we need them and don't want to set them in the conventional way
                    for a in state_attrs:
                        vals.pop(attrs.index(a))
                        attrs.pop(attrs.index(a))
                    
                    #: The string representation of the index (1-based)
                    I = str(i+1)
                    
                    #Temperature and pressure provided
                    if 'injection_state_temp_'+I in state_attrs and 'injection_state_pressure_'+I in state_attrs:
                        injection_temp = state_vals[state_attrs.index('injection_state_temp_'+I)]
                        injection_pressure = state_vals[state_attrs.index('injection_state_pressure_'+I)]
                        self.Lines[i].state.set_state(inletState.Fluid,
                                                      T=injection_temp, 
                                                      P=injection_pressure)
                        
                    #Dew temperature and superheat provided
                    elif 'injection_state_sat_temp_'+I in state_attrs and 'injection_state_superheat_'+I in state_attrs:
                        injection_sat_temp = state_vals[state_attrs.index('injection_state_sat_temp_'+I)]
                        injection_superheat = state_vals[state_attrs.index('injection_state_superheat_'+I)]
                        injection_temp = injection_sat_temp + injection_superheat
                        import CoolProp.CoolProp as CP
                        injection_pressure = CP.PropsSI('P','T',injection_sat_temp,'Q',1.0,inletState.Fluid)/1000.0
                        self.Lines[i].state.set_state(inletState.Fluid,
                                                      T=injection_temp, 
                                                      P=injection_pressure)
                        
                    else:
                        raise ValueError('Invalid combination of injection states: '+str(state_attrs))
                    
                elif num_inj_state_params == 1:
                    import textwrap
                    string = textwrap.dedent(
                             """
                             Sorry but you need to provide two variables for the injection
                             state in parametric table to fix the state.  
                             
                             If you want to just modify the saturated temperature, add the superheat as a
                             variable and give it one element in the parametric table
                             """
                             )
                    dlg = wx.MessageDialog(None,string)
                    dlg.ShowModal()
                    dlg.Destroy()
                    raise ValueError('Must provide two state variables in the parametric table for injection line')
                    
                elif num_inj_state_params >2:
                    raise ValueError ('Only two inlet state parameters can be provided in parametric table')
            
            return attrs,vals
        
        
        def apply_port_terms(attrs,vals):
            phi_params = [(par,val) for par, val in zip(attrs,vals) if par.startswith('injection_phi')]
            num_phi_params = len(phi_params)        
            
            if num_phi_params > 0:
                #Unzip the parameters (List of tuples -> tuple of lists)
                phi_attrs, phi_vals = zip(*phi_params)
                
                # Remove all the entries that correspond to the angles 
                # we need them and don't want to set them in the conventional way
                for a in phi_attrs:
                    i = attrs.index(a)
                    vals.pop(i)
                    attrs.pop(i)
                    
                for attr,val in zip(phi_attrs, phi_vals):
    
                    # Term might look like something like 'injection_phi_1_2'
                    # i would be 0, j would be 1
                    #indices are zero-based
                    j = int(attr.rsplit('_',1)[1])-1
                    i = int(attr.rsplit('_',2)[1])-1
                    
                    self.Lines[i].ports_list[j].phi_inj_port.SetValue(str(val))
            
            return attrs,vals
        
        #Apply all the line terms and get back the lists
        attrs, vals = apply_line_terms(attrs,vals)
        #Apply all the line terms and get back the lists
        attrs, vals = apply_port_terms(attrs,vals)
    
        return attrs,vals
      

VI_template = """           
        
from PDSim.scroll.core import Port
import scipy.interpolate

phi_list = {phi_list:s}
involute_list = {involute_list:s}
offset_list = {offset_list:s}
D_list = {D_list:s}
X_d_list = {X_d_list:s}
X_d_backflow_list = {X_d_backflow_list:s}

port_tube_list = {port_tube_list:s}

Ltube_list = {Ltube_list:s}
IDtube_list = {IDtube_list:s}
Ttube_list = {Ttube_list:s}
RHOtube_list = {RHOtube_list:s}

sim.fixed_scroll_ports = []

for phi,involute,offset,D,parent,X_d,X_d_backflow in zip(phi_list,involute_list,offset_list,D_list,port_tube_list,X_d_list,X_d_backflow_list):
    p = Port()
    p.phi = phi
    p.involute = involute
    p.offset = offset
    p.D = D
    p.parent = parent
    p.X_d = X_d
    p.X_d_backflow = X_d_backflow
    sim.fixed_scroll_ports.append(p)
    
print mdot_guess*0.1
for i,(L,ID,T,D) in enumerate(zip(Ltube_list,IDtube_list,Ttube_list,RHOtube_list)):
    sim.add_tube(Tube(key1 = 'VITube'+str(i+1)+'.1',
                      key2 = 'VITube'+str(i+1)+'.2',
                      L = L,
                      ID = ID,
                      mdot = mdot_guess*0.1, 
                      State1 = State.State(inletState.Fluid,
                                           dict(T = T, D = D)
                                           ),
                      fixed = 1,
                      TubeFcn = sim.TubeCode
                      )     
                  )
        
##This block can be uncommented to plot the scroll wrap with the VI ports
##plotScrollSet(0, geo = sim.geo)
#for port in sim.VI_ports:
#    x,y = scroll_geo.coords_inv(port.phi, sim.geo, 0, 'f'+port.involute)
#    nx,ny = scroll_geo.coords_norm(port.phi, sim.geo, 0, 'f'+port.involute)
#    
#    x0 = x - nx*port.offset
#    y0 = y - ny*port.offset
#    
#    t = np.linspace(0, 2*pi)
#    plt.plot(port.D/2.0*np.cos(t) + x0, port.D/2.0*np.sin(t) + y0,'b')
#    
#plt.savefig('wrap_with_VI_ports.png',transparent=True)
#plt.close()
#quit()
        
#  Calculate the areas between each VI port and every control volume
sim.calculate_port_areas()
    
for port in sim.fixed_scroll_ports:
    for partner in port.area_dict:
    
        #  Create a spline interpolator object for the area between port and the partner chamber
        A_interpolator = scipy.interpolate.splrep(port.theta, port.area_dict[partner], k = 2, s = 0)
        
        #  Add the flow between the EVI control volume and the chamber through the port
        sim.add_flow(FlowPath(key1 = 'VITube'+str(port.parent+1)+'.2',
                              key2 = partner,
                              MdotFcn = sim.INTERPOLATING_NOZZLE_FLOW,
                              MdotFcn_kwargs = dict(X_d = port.X_d,
                                                    X_d_backflow = port.X_d_backflow,
                                                    upstream_key = 'VITube'+str(port.parent+1)+'.2',
                                                    A_interpolator = A_interpolator
                                                    )
                             )
                     )
"""
class ScrollInjectionPlugin(pdsim_plugins.PDSimPlugin):
    """
    A plugin that adds the injection ports for the scroll compressor
    """
    
    #: A short description of the plugin 
    short_description = 'Vapor injection for scroll compressor'
        
    def should_enable(self):
        """
        Only enable if it is a scroll type compressor
        """
        if self.GUI.machinefamily.lower().startswith('scroll'):
            return True
        else:
            return False
    
    def get_config_chunk(self):
        """
        The chunk for the configuration file
        """
        
        chunk = []
        for line in self.injection_panel.Lines:
            l = {}
            l['Length'] = float(line.Lval.GetValue())
            l['ID'] = float(line.IDval.GetValue())
            State = line.state.GetState()
            l['inletState'] = dict(T = State.T, rho = State.rho, Fluid = State.Fluid)
            
            #Then the things related to the ports for this line
            l['ports'] = []
            for portpanel in line.ports_list:

                port = {}
                
                #The things that are for the port
                port['phi'] = float(portpanel.phi_inj_port.GetValue())
                inv = portpanel.involute.GetStringSelection()
                if inv == 'Inner involute':
                    port['involute'] = 'i'
                elif inv == 'Outer involute':
                    port['involute'] = 'o'
                else:
                    raise ValueError
                
                port['check_valve'] = portpanel.check_valve.IsChecked()
                port['symmetric'] = portpanel.SymmTarget.GetStringSelection()
                port['D'] = float(portpanel.D_port.GetValue())
                port['offset'] = float(portpanel.offset_port.GetValue())
                port['X_d'] = float(portpanel.X_d.GetValue())
                port['X_d_backflow'] = float(portpanel.X_d_backflow.GetValue())
                
                l['ports'].append(port)
            
            chunk.append(l)
        
        top_level_dict = {'Plugin:ScrollInjectionPlugin':chunk}
        return top_level_dict    
        
    def activate(self, event = None, config = ''):
        """
        Called to activate or deactivate the injection plugin
        """
        #: The inputs toolbook that contains all the input panels
        ITB = self.GUI.MTB.InputsTB
        
        if not self._activated:
            #Add the panel to the inputs panel
            #name is the internal name, also used in saving and loading 
            # config files
            self.injection_panel = InjectionInputsPanel(ITB, name = 'Plugin:ScrollInjectionPlugin')
            ITB.AddPage(self.injection_panel,"Injection")
            self._activated = True
            self.injection_panel.build_from_configfile(config)
        else:
            page_names = [ITB.GetPageText(I) for I in range(ITB.GetPageCount())]
            I = page_names.index("Injection")
            ITB.RemovePage(I)
            self.injection_panel.Destroy()
            del self.injection_panel
            self._activated = False
        
    def get_script_chunks(self):
        """
        Add the necessary things for the scroll compressor injection
        """
        
        Ltube_list,IDtube_list,Ttube_list,RHOtube_list = [],[],[],[]
        phi_list,involute_list,offset_list,check_valve_list, D_list, port_tube_list,X_d_list,X_d_backflow_list = [],[],[],[],[],[],[],[]
            
        #IEPs are children of injection_panel that are instances of InjectionElementPanel class
        IEPs = [child for child in self.injection_panel.scrolled_panel.Children if isinstance(child,InjectionElementPanel)]
        for i, IEP in enumerate(IEPs):
            Ltube_list.append(float(IEP.Lval.GetValue()))
            IDtube_list.append(float(IEP.IDval.GetValue()))
            injState = IEP.state.GetState()
            Ttube_list.append(injState.T)
            RHOtube_list.append(injState.rho)
            
            Ports = [c for c in IEP.Children if isinstance(c,InjectionPortPanel)]
            for j,child in enumerate(Ports):
                v = child.get_values()
                phi_list.append(v['phi'])
                involute_list.append(v['involute'])
                check_valve_list.append(v['check_valve'])
                D_list.append(v['D'])                
                offset_list.append(v['offset'])
                port_tube_list.append(i)
                X_d_list.append(v['X_d'])
                X_d_backflow_list.append(v['X_d_backflow'])
                
                                        
        post_build = VI_template.format(**{k:str(v) for k,v in locals().iteritems()})
        return dict(post_build = post_build)
                        
#    def post_process(self, sim):
#        """
#        Post-process the results from the simulation in order to calculate any parameters that
#        are required
#        
#        This function will be called by OnIdle in GUI Main frame when run finishes
#        """
#
#        sim.injection.massflow={}
#        #:the ratio of the injection flow rate to the suction flow rate
#        sim.injection.flow_ratio={}
#        #injection pressure
#        sim.injection.pressure={}
#        #injection temperature
#        sim.injection.temperature={}
#        
#        #The tubes that are injection tubes have a key1 that starts with 'injection_line'
#        ITubes = [T for T in sim.Tubes if T.key1.startswith('injection_line')]
#        
#        for i,Tube in enumerate(ITubes):
#            key = Tube.key1
#            sim.injection.massflow[i+1]=sim.FlowsProcessed.mean_mdot[key]
#            sim.injection.flow_ratio[i+1]=(sim.injection.massflow[i+1]/
#                                           sim.mdot)
#            sim.injection.pressure[i+1] = Tube.State1.p
#            sim.injection.temperature[i+1] = Tube.State1.T
#                
#        #Save a local copy of a pointer to the simulation
#        self.simulation = sim

if __name__=='__main__':
    print 'I am not meant to be run directly'
        
        
        