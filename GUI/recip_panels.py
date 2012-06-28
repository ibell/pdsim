import pdsim_panels
import wx
from math import pi

class GeometryPanel(pdsim_panels.PDPanel):
    """
    The geometry panel of the reciprocating compressor
    Loads all parameters from the configuration file
    """
    def __init__(self,parent,configfile,**kwargs):
        pdsim_panels.PDPanel.__init__(self,parent,**kwargs)
        
        #Loads all the parameters from the config file
        configdict,descdict = self.get_from_configfile(configfile,'GeometryPanel')
        
        # Things in self.items are linked through to the module code where 
        # it attempts to set the attribute.  They are also automatically
        # written to configuration file
        self.items = [
        dict(attr='piston_diameter'),
        dict(attr='piston_length'),
        dict(attr='crank_length'),
        dict(attr='connecting_rod_length'),
        dict(attr='x_TDC', tooltip='The distance from the top of the cylinder to the piston head at top dead center'),
        dict(attr='V_backchamber'),
        ]
        
        sizer = wx.FlexGridSizer(cols=2, vgap=4, hgap=4)
        
        self.ConstructItems(self.items,sizer,configdict,descdict)
            
        self.SetSizer(sizer)
        sizer.Layout()
        
class MassFlowPanel(pdsim_panels.PDPanel):
    
    def __init__(self, parent, configfile,**kwargs):
    
        pdsim_panels.PDPanel.__init__(self, parent,**kwargs)
        
        #Loads all the parameters from the config file
        self.configdict,self.descdict = self.get_from_configfile(configfile,'MassFlowPanel')
        
        self.items1 = [
        dict(attr='d_discharge'),
        dict(attr='d_suction'),
        ]
        box_sizer = wx.BoxSizer(wx.VERTICAL)
        box_sizer.Add(wx.StaticText(self,-1,"Required Inputs"))
        box_sizer.Add(wx.StaticLine(self,-1,(25, 50), (300,1)))
        
        sizer = wx.FlexGridSizer(cols=2, vgap=4, hgap=4)
        self.ConstructItems(self.items1,sizer,self.configdict,self.descdict)

        box_sizer.Add(sizer)
        box_sizer.Add((10,10))
        box_sizer.Add(wx.StaticText(self,-1,"Valve Inputs"))
        box_sizer.Add(wx.StaticLine(self,-1,(25, 50), (300,1)))
        
        self.items2 = [
        dict(attr='valve_E'),
        dict(attr='valve_d'),
        dict(attr='valve_h'),
        dict(attr='valve_l'),
        dict(attr='valve_a'),
        dict(attr='valve_x_stopper'),
        dict(attr='valve_rho'),
        dict(attr='valve_C_D'),
        ]
        sizer = wx.FlexGridSizer(cols=2, vgap=4, hgap=4)
        self.ConstructItems(self.items2,sizer,self.configdict,self.descdict)
        box_sizer.Add(sizer)
        self.SetSizer(box_sizer)
        box_sizer.Layout()
        
        self.items=self.items1+self.items2

        
class MechanicalLossesPanel(pdsim_panels.PDPanel):
    
    def __init__(self, parent, configfile,**kwargs):
    
        pdsim_panels.PDPanel.__init__(self, parent,**kwargs)
        
        #Loads all the parameters from the config file (case-sensitive)
        self.configdict,self.descdict = self.get_from_configfile(configfile,'MechanicalLossesPanel')
        
        self.items = [
        dict(attr='eta_motor'),
        dict(attr='h_shell'),
        dict(attr='A_shell'),
        dict(attr='Tamb'),
        dict(attr='mu_oil'),
        dict(attr='delta_gap'),
        ]
        
        sizer = wx.FlexGridSizer(cols=2, vgap=4, hgap=4)
        
        self.ConstructItems(self.items,sizer,self.configdict,self.descdict)

        self.SetSizer(sizer)
        sizer.Layout()
        
class StatePanel(pdsim_panels.PDPanel):
    
    def __init__(self, parent, configfile,**kwargs):
    
        pdsim_panels.PDPanel.__init__(self, parent,**kwargs)
        
        #Loads all the parameters from the config file (case-sensitive)
        self.configdict, self.descdict = self.get_from_configfile(configfile,'StatePanel')
        
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
        self.SuctionState = pdsim_panels.StatePanel(self,Fluid=Fluid,T=T,rho=rho)
        box_sizer.Add(self.SuctionState)
        
        box_sizer.Add((20,20))
        box_sizer.Add(wx.StaticText(self,-1,"Discharge State"))
        box_sizer.Add(wx.StaticLine(self,-1,(25, 50), (300,1)))
        
        self.cmbDischarge = wx.ComboBox(self)
        self.cmbDischarge.AppendItems(['Discharge pressure [kPa]', 'Pressure ratio [-]'])
        self.cmbDischarge.SetStringSelection(self.Discharge_key)
        self.DischargeValue = wx.TextCtrl(self, value = self.Discharge_value)
        self.cmbDischarge.Bind(wx.EVT_COMBOBOX, self.OnChangeDischarge)
        
        sizer = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
        sizer.AddMany([self.cmbDischarge, self.DischargeValue])
        box_sizer.Add(sizer)
        
        self.SetSizer(box_sizer)
        sizer.Layout()
        
    def OnChangeDischarge(self, event):
        p_suction = self.SuctionState.GetState().p

        if self.cmbDischarge.GetStringSelection() == 'Discharge pressure [kPa]':
            pratio = float(self.DischargeValue.GetValue())
            p = pratio*p_suction
            self.DischargeValue.SetValue(str(p))
        elif self.cmbDischarge.GetStringSelection() == 'Pressure ratio [-]':
            p_disc = float(self.DischargeValue.GetValue())
            pratio = p_disc/p_suction
            self.DischargeValue.SetValue(str(pratio))
        else:
            raise KeyError
        
    def post_get_from_configfile(self, key, value):
        Dummy, value, key = value.split(',')
        self.Discharge_key = key
        self.Discharge_value = str(value)
        
    def post_calculate(self, recip):
        recip.inletState = self.SuctionState.GetState()
        if self.cmbDischarge.GetStringSelection() == 'Discharge pressure [kPa]':
            recip.discharge_pressure = float(self.DischargeValue.GetValue())
        elif self.cmbDischarge.GetStringSelection() == 'Pressure ratio [-]':
            p_suction = self.SuctionState.GetState().p
            p_ratio = float(self.DischargeValue.GetValue())
            recip.discharge_pressure = p_ratio * p_suction
        
    def post_prep_for_configfile(self):
        """
        Write a string representation of the state
        """
        State_ = self.SuctionState.GetState()
        StateString = 'inletState = State,'+State_.Fluid+','+str(State_.T)+','+str(State_.rho)
        DischargeString = 'discharge = Discharge,'+str(self.DischargeValue.GetValue())+','+self.cmbDischarge.GetStringSelection()
        return StateString+'\n'+DischargeString