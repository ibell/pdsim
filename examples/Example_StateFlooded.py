from PDSim.core.state_flooded import StateFlooded

S = StateFlooded('Nitrogen', 'POE', 101.325, 300, 0.5, 'HEM')
S.update(dict(T=300, P=101.325, xL=0.0))
print 'Nitrogen'
print S
print '--------------------------'
S = StateFlooded('R245FA', 'POE', 101.325, 300, 0.5, 'HEM')
S.update(dict(T=300, P=101.325, xL=0.3))
print 'R245FA'
print S
print '--------------------------'

from CoolProp.CoolProp import State
S = State('Nitrogen', dict(T=300, P=101.325, xL=0.))
print 'CoolProp StateClass'
print S
print '--------------------------'