import scipy
from ReactionNetworks import *
#WARNING: This must match the synthetic experiment name
net = Network('ActivationInhibition')

net.add_compartment('cell')

net.add_parameter('k0', 1, name = r'k_{0}')
net.add_parameter('k1', 1, name = r'k_{1}')
net.add_parameter('k2', 1, name = r'k_{2}')
net.add_parameter('k3', 1, name = r'k_{3}')
net.add_parameter('k5',  1, name = r'k5')
net.add_parameter('k6',  1, name = r'k6')
net.add_parameter('ET', 1, name = r'E_T')
net.add_parameter('km5',.1, name = r'k_{-3}')
net.add_parameter('km6',.1, name = r'k_{-4}')
net.add_parameter('k4', .5, name = r'k_{4}')
net.add_parameter('S', .5, name = r'k_{4}')



net.add_species('X', 'cell', 1)
net.add_species('R', 'cell', 1)

net.add_parameter('u1')
net.add_assignment_rule('u1', 'k5*R')
net.add_parameter('u2')
net.add_assignment_rule('u2', 'k6')
net.add_parameter('J1')
net.add_assignment_rule('J1', 'km5/ET')
net.add_parameter('J2')
net.add_assignment_rule('J2', 'km6/ET')

net.add_parameter('B')
net.add_assignment_rule('B', 'u2 - u1 + u2*J1 + u1*J2')

net.add_parameter('G')
net.add_assignment_rule('G', '2*u1*J2/ (B + sqrt(B**2 - 4*(u2-u1)*u1*J2))')
net.add_parameter('Ep')
net.add_assignment_rule('Ep', 'ET*G')

net.add_rate_rule('R', 'k0*Ep + k1*S - k2*X*R')
net.add_rate_rule('X', 'k3*Ep - k4*X ')


network = net
int_time = (0, 1000)
