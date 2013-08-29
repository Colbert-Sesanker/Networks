import scipy
from ReactionNetworks import *
# WARNING: This must match the synthetic experiment name
net = Network('goodwin')

net.add_compartment('cell')

net.add_parameter('k1', 1)
net.add_parameter('k2', 1)
net.add_parameter('k3', 1)
net.add_parameter('k4', .2)
net.add_parameter('k5', .15)
net.add_parameter('k6', .1)
net.add_parameter('k', 1)
net.add_parameter('n', 9)

net.add_species('x', 'cell', .0348)
net.add_species('y', 'cell', .347)
net.add_species('z', 'cell', 1.735)

net.add_rate_rule('x', 'k1 / (k + z**n) - k4*x')
net.add_rate_rule('y', 'k2*x - k5*y')
net.add_rate_rule('z', 'k3*y - k6*z')

network = net
int_time = (0, 600)
