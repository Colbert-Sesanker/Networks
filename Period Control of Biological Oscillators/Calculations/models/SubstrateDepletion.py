import scipy
from ReactionNetworks import *
#WARNING: This must match the synthetic experiment name
net = Network('SubstrateDepletion')

net.add_compartment('cell')

net.add_parameter('k0', 1, name = r'k_{0}')
net.add_parameter('k1', 1, name = r'k_{1}')
net.add_parameter('k2', 1, name = r'k_{2}')
net.add_parameter('k3', 1, name = r'k_{3}')
net.add_parameter('k4', 1, name = r'k_{4}')
net.add_parameter('S',  1, name = r'S')
net.add_parameter('ET', 1, name = r'E_T')
net.add_parameter('km3',.1, name = r'k_{-3}')
net.add_parameter('km4',.1, name = r'k_{-4}')
net.add_parameter('k0p',.1, name = r'k_{0p}')



net.add_species('X', 'cell', 1)
net.add_species('R', 'cell', 1)

net.add_parameter('u1')
net.add_assignment_rule('u1', 'k3*R')
net.add_parameter('u2')
net.add_assignment_rule('u2', 'k4')
net.add_parameter('J1')
net.add_assignment_rule('J1', 'km3/ET')
net.add_parameter('J2')
net.add_assignment_rule('J2', 'km4/ET')

net.add_parameter('B')
net.add_assignment_rule('B', 'u2 - u1 + u2*J1 + u1*J2')

net.add_parameter('G')
net.add_assignment_rule('G', '2*u1*J2/ (B + sqrt(B**2 - 4*(u2-u1)*u1*J2))')
net.add_parameter('Ep')
net.add_assignment_rule('Ep', 'ET*G')

net.add_rate_rule('X', 'k1*S - (k0p + k0*Ep)*X')
net.add_rate_rule('R', '(k0p + k0*Ep)*X - k2*R')

network = net
int_time = (0, 100)
