import os 

from PDSim.scroll.core import Scroll
from scroll_compressor import Compressor

from CoolProp import CoolProp as CP
CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, os.getenv('RPPREFIX'))

# CP.set_debug_level(1)

# Use this one to use R410A as a pseudo-pure fluid
Compressor(Scroll, Ref='R410A', HDF5file='R410APPF.h5')

# Use this one to use R410A as a mixture w/ REFPROP
Compressor(Scroll, Ref='REFPROP::R32[0.697614699375863]&R125[0.302385300624138]', HDF5file='REFPROPR410Amix.h5')

# Use this one to use R410A as a mixture w/ REFPROP and bicubic interpolation
Compressor(Scroll, Ref='BICUBIC&REFPROP::R32[0.697614699375863]&R125[0.302385300624138]', HDF5file='REFPROPBICUBICR410Amix.h5')

# Use this one to use R410A as a mixture w/ the HEOS backend of CoolProp
Compressor(Scroll, Ref='HEOS::R32[0.697614699375863]&R125[0.302385300624138]', HDF5file='HEOSR410Amix.h5')