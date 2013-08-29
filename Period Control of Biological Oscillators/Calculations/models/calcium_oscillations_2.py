import scipy
from ReactionNetworks import *
# WARNING: This must match the synthetic experiment name
net = Network('calcium_oscillations_2')

net.add_compartment('cell')

net.add_parameter('c1', 6.64)
net.add_parameter('c2', 5)
net.add_parameter('c3', .0000313)
net.add_parameter('c4',  .435)
net.add_parameter('c5',  2)
net.add_parameter('c6', .5)
net.add_parameter('c7', .6)
net.add_parameter('c8',  6.64)
net.add_parameter('c9', 5)
net.add_parameter('c10',.0000313)
net.add_parameter('k1', .1)
net.add_parameter('k2', .15)
net.add_parameter('k3', 1)
net.add_parameter('k4', .1)
net.add_parameter('k5', .15)
net.add_parameter('n', 3)
net.add_parameter('m', 2)
net.add_parameter('k', 2)
net.add_parameter('p', 3.3)


net.add_species('x', 'cell', .1)
net.add_species('y', 'cell', .01)
net.add_species('z', 'cell', 25)

net.add_rate_rule('x', '  (z * c1 * y**3) / (k1 + y)**n   -   (c2 * x**2) / (x + k2)**m + c3*z**k - c6*(x / c7)**(p) + c6')
net.add_rate_rule('y', '  (c4 * x) / (k3 + x)  - c5*y')
net.add_rate_rule('z', '- (z * c8 * y**3) / (k4 + y)**n   +  (c9 * x**2) / (x + k5)**m - c10*z**k')
network = net
int_time = (0, 300)
