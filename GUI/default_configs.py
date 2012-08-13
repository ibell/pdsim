# -*- coding: latin-1 -*-

import StringIO

recip_template = (
"""
[Globals]
Type = recip
Mode = compressor

[GeometryPanel]
piston_diameter = float,0.02,Piston diameter [m]
piston_length = float,0.02,Piston length [m]
crank_length = float,0.01,Crank length [m]
connecting_rod_length = float,0.04,Connecting rod length [m]
dead_volume_perc = float,4.0,Dead volume percentage [%%]
x_TDC = float,0.005,Distance to piston at TDC [m]
shell_volume = float,100e-6,Shell volume [m³]

[MassFlowPanel]
d_discharge = float,0.0059,Discharge port diameter [m]
d_suction = float,0.0059,Suction port diameter [m]
valve_E = float,1.93e+11,Youngs Modulus [Pa]
valve_d = float,0.007,Valve diameter [m]
valve_h = float,0.0001532,Valve thickness [m]
valve_l = float,0.018,Valve length [m]
valve_a = float,0.014,Valve distance from anchor [m]
valve_x_stopper = float,0.0018,Valve distance to stopper [m]
valve_rho = float,8000.0,Valve metal density [kg/m³]
valve_C_D = float,1.17,Valve drag coefficient [-]

[MechanicalLossesPanel]
eta_motor = float,0.95,Motor efficiency [-]
h_shell = float,0.01,Shell air-side heat transfer coefficient [kW/m²/K]
A_shell = float,0.040536,Shell Area [m²]
Tamb = float,298.0,Ambient temperature [K]
mu_oil = float,0.0086,Oil viscosity [Pa-s]
delta_gap = float,2e-05,Gap width [m]

[StatePanel]
omega = float,377.0,Rotational speed [rad/s]
inletState = State,R404A,283.15,5.75
discharge = Discharge,2.0,Pressure ratio [-]

[ParametricPanel]
Term1 = Term,Rotational speed [rad/s],250.0;300.0

[SolverInputsPanel]
Cycle = Cycle,Euler,7000

[OutputDataPanel]
selected  = selected,['run_index'; 'mdot'; 'eta_v'; 'eta_oi'; 'Td']
"""
)

scroll_template=(
"""
[Globals]
Type = scroll
Mode = compressor

[GeometryPanel]
Vdisp = float,104.8e-6,Displacement volume / revolution [m³/rev]
Vratio = float,1.61,Built-in volume ratio [-]
ro = float, 0.005, Orbiting radius [m]
t = float,0.004,Scroll wrap thickness [m]
delta_flank = float,3e-6,Flank gap width [m]
delta_radial = float,3e-6,Radial gap width [m]

[MassFlowPanel]
d_discharge = float,0.0059,Discharge port diameter [m]
d_suction = float,0.0059,Suction port diameter [m]
inlet_tube_length = float,0.3,Inlet tube length [m]
inlet_tube_ID = float,0.02,Inlet tube inner diameter [m]
outlet_tube_length = float,0.3,Outlet tube length [m]
outlet_tube_ID = float, 0.02, Outlet tube inner diameter [m]

[MechanicalLossesPanel]
eta_motor = float,0.95,Motor efficiency [-]
h_shell = float,0.01,Shell air-side heat transfer coefficient [kW/m²/K]
A_shell = float,0.040536,Shell Area [m²]
Tamb = float,298.0,Ambient temperature [K]
mu_oil = float,0.0086,Oil viscosity [Pa-s]
delta_gap = float,2e-05,Gap width [m]

[StatePanel]
omega = float,377.0,Rotational speed [rad/s]
inletState = State,R410A,283.15,5.75
discharge = Discharge,2.0,Pressure ratio [-]

[ParametricPanel]
Term1 = Term,Rotational speed [rad/s],250.0;300.0

[SolverInputsPanel]
Cycle = Cycle,Heun,7000

[OutputDataPanel]
selected  = selected,['run_index'; 'mdot'; 'eta_v'; 'eta_a']
"""
)

def get_recip_defaults():
    """
    Create a cStringIO object that can be read by the configparser
    """
    return StringIO.StringIO(recip_template)

def get_scroll_defaults():
    """
    Create a cStringIO object that can be read by the configparser
    """
    return StringIO.StringIO(scroll_template)
    