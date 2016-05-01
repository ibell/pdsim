from PDSim.core.state_flooded import StateFlooded

S=StateFlooded("R245fa", 'POE', 1200, 373.15, 0.1, 'HEM')
print S

S = StateFlooded('Nitrogen', 'POE', 101.325, 300, 0.5, 'HEM')
S.update(dict(T=300, P=101.325, xL=1))
print 'Nitrogen'
print S.get_h()
print '--------------------------'
S = StateFlooded('R245FA', 'POE', 101.325, 300, 0.5, 'HEM')
S.update(dict(T=300, P=101.325, xL=0.3))
print 'R245FA'
print S.get_kstar()
print '--------------------------'

from CoolProp.CoolProp import State
S = State('Nitrogen', dict(T=300, P=101.325, xL=0.))
print 'CoolProp StateClass'
print S
print '--------------------------'