from scroll_compressor import Compressor
import CoolProp.CoolProp as CP
#CP.set_debug_level(1)

# Use this one to use R410A as a mixture
Compressor(Ref='BICUBIC&REFPROP::R32[0.697614699375863]&R125[0.302385300624138]')

# Use this one to use R410A as a pure fluid
#Compressor(Ref='R410A')