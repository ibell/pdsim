from __future__ import print_function

from PDSim.scroll.core import Scroll
from scroll_compressor import Compressor

# from CoolProp import CoolProp as CP
# CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, r'/home/ian/REFPROP911')

# CP.set_debug_level(1)

# Use this one to use R410A as a mixture
Compressor(Scroll, Ref='BICUBIC&REFPROP::R32[0.697614699375863]&R125[0.302385300624138]')

# Use this one to use R410A as a pure fluid
#Compressor(Ref='R410A')
