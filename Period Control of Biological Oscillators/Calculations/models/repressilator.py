import scipy
from ReactionNetworks import *
# WARNING: This must match the synthetic experiment name
net = Network('repressilator')

net.add_compartment('cell')

net.add_parameter('k1', 2.35804)
net.add_parameter('k2', 4.42269)
net.add_parameter('k3', 4.80922)
net.add_parameter('n1', 5.03005)
net.add_parameter('n2', 5.73448)
net.add_parameter('n3', 6.05167)
net.add_parameter('a1', 5.73702)
net.add_parameter('a2', 6.92108)
net.add_parameter('a3', 7.46407)
net.add_parameter('b1', .3284)
net.add_parameter('b2', .4967)
net.add_parameter('b3', .4518)
net.add_parameter('y1', .0908)
net.add_parameter('y2', .08093)
net.add_parameter('y3', .11444)



net.add_species('A', 'cell', .48199)
net.add_species('B', 'cell', 5.11385 )
net.add_species('C', 'cell', 105.77422)


net.add_func_def('hill_promote', ('x', 'k', 'n'), 'x**n / (k**n + x**n)')
net.add_func_def('hill_repress', ('x', 'k', 'n'), 'k**n / (k**n + x**n)')


net.add_rate_rule('A', 'a1*hill_repress(B, k1, n1) - b1*A + y1')
net.add_rate_rule('B', 'a2*hill_repress(C, k2, n2) - b2*B + y2')
net.add_rate_rule('C', 'a3*hill_repress(A, k3, n3) - b3*C + y3')

network = net
int_time = (0, 200)
