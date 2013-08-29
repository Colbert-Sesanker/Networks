import scipy
from ReactionNetworks import *
# WARNING: This must match the synthetic experiment name
net = Network('calcium_oscillations')

net.add_compartment('cell')

net.add_parameter('c1', 6.64)
net.add_parameter('c2', 5)
net.add_parameter('c3', .0004)#.0000313)
net.add_parameter('c4', .435)
net.add_parameter('c5',  2)
net.add_parameter('c6', .5)
net.add_parameter('c7', .6)
net.add_parameter('k1', .1)
net.add_parameter('k2', .15)
net.add_parameter('k3', 1)


net.add_species('x', 'cell', .1)
net.add_species('y', 'cell', .01)
net.add_species('z', 'cell', 25)

net.add_rate_rule('x', '  (z * c1 * y**3) / (k1 + y)**3   -  (c2 * x**2) / (x + k2)**2 + c3*z**2 - c6*(x / c7)**(3.3) + c6')
net.add_rate_rule('y', '  (c4 * x) / (k3 + x)  - c5*y')
net.add_rate_rule('z', '- (z * c1 * y**3) / (k1 + y)**3   +  (c2 * x**2) / (x + k2)**2 - c3*z**2')
network = net
int_time = (0, 300)
