
function dy = Repressilator_Positive_Sens(t, y, p)

% Initial conditions for FHN model are V = -1, R = 1

% Set up species
A   = y(1);
B   = y(2);
C   = y(3);
Aa  = y(1);
Ba  = y(2);
Ca  = y(3);
An  = y(1);
Bn  = y(2);
Cn  = y(3);
Ak  = y(1);
Bk  = y(2);
Ck  = y(3);
Ra = y(4);
Vb = y(5);
Rb = y(6);
Vc = y(7);
Rc = y(8);
Vaa = y(9);
Raa = y(10);
Vab = y(11);
Rab = y(12);
Vac = y(13);
Rac = y(14);
Vbb = y(15);
Rbb = y(16);
Vbc = y(17);
Rbc = y(18);
Vcc = y(19);
Rcc = y(20);

% Set up parameters
k1 = p(1);
k2 = p(2);
k3 = p(3);
k4 = p(4);
n1 = p(5);
n2 = p(6);
n3 = p(7);
a1 = p(8);
a2 = p(9);
a3 = p(10);
a4 = p(11);
b1 = p(12);
b2 = p(13);
b3 = p(14);
y1 = p(15);
y2 = p(16);
y3 = p(17);

% Define hill functions
hill_promote = @(x, k, n)  x^n / (k^n + x^n);
hill_repress = @(x, k, n)  k^n / (k^n + x^n);

% Evaluate equations
dy    = zeros(2,1);    % a column vector

dy(1) = a4*hill_promote(A, k4, n4) + a1*hill_repress(B, k1, n1) - b1*A + y1;
dy(2) = a2*hill_repress(C, k2, n2) - b2*B + y2;
dy(3) = a3*hill_repress(B, k3, n3) - b2*C + y3;

dy(2) = -(V-a+b*R)/c;

dy(3) = (c-c*V^2)*Va + c*Ra;
dy(4) = (-1/c)*Va + (-b/c)*Ra + 1/c;
dy(5) = (c-c*V^2)*Vb + c*Rb;
dy(6) = (-1/c)*Vb + (-b/c)*Rb - R/c;
dy(7) = (c-c*V^2)*Vc + c*Rc + V - (V^3)/3 + R;
dy(8) = (-1/c)*Vc + (-b/c)*Rc + (V-a+b*R)/(c^2);

dy(9) = (-2)*V*c*Va*Va - c*(V^2 - 1)*Vaa + c*Raa;
dy(10) = -1/c*Vaa - b/c*Raa;
dy(11) = (-2)*V*c*Vb*Va - c*(V^2 - 1)*Vab + c*Rab;
dy(12) = -1/c*Vab - 1/c*Ra - b/c*Rab;
dy(13) = (-2)*V*c*Vc*Va + 1 - V^2*Va - c*(V^2 - 1)*Vac + 1*Ra + c*Rac;
dy(14) = 1/c^2*Va - 1/c*Vac + b/c^2*Ra - b/c*Rac - 1/c^2;
dy(15) = (-2)*V*c*Vb*Vb - c*(V^2 - 1)*Vbb + c*Rbb;
dy(16) = -1/c*Vbb - 1/c*Rb - b/c*Rbb - 1/c*Rb;
dy(17) = (-2)*V*c*Vc*Vb + 1 - V^2*Vb - c*(V^2 - 1)*Vbc + 1*Rb + c*Rbc;
dy(18) = 1/c^2*Vb - 1/c*Vbc + b/c^2*Rb - b/c*Rbc - 1/c*Rc + R/c^2;
dy(19) = (-2)*V*c*Vc*Vc + 1 - V^2*Vc - c*(V^2 - 1)*Vcc + 1 - V^2*Vc + 1*Rc + c*Rcc + 1*Rc;
dy(20) = 1/c^2*Vc - 1/c*Vcc + 1/c^2*Vc + b/c^2*Rc - b/c*Rcc + b/c^2*Rc - (2*(V - a + R*b))/c^3;


end




net.add_parameter('k1', 2.35804)
net.add_parameter('k2', 4.42269)
net.add_parameter('k3', 4.80922)
net.add_parameter('k4', 5)
net.add_parameter('n1', 5.03005)
net.add_parameter('n2', 5.73448)
net.add_parameter('n3', 6.05167)
net.add_parameter('n4', 7)
net.add_parameter('a1', 5.73702)
net.add_parameter('a2', 6.92108)
net.add_parameter('a3', 7.46407)
net.add_parameter('a4', .7)
net.add_parameter('b1', .3284)
net.add_parameter('b2', .4967)
net.add_parameter('b3', .4518)
net.add_parameter('y1', .908)
net.add_parameter('y2', .8093)
net.add_parameter('y3', 1.1444)



net.add_species('A', 'cell', .48199)
net.add_species('B', 'cell', 5.11385 )
net.add_species('C', 'cell', 105.77422)


net.add_func_def('hill_promote', ('x', 'k', 'n'), 'x**n / (k**n + x**n)')
net.add_func_def('hill_repress', ('x', 'k', 'n'), 'k**n / (k**n + x**n)')


net.add_rate_rule('A', 'a1*hill_repress(B, k1, n1) - b1*A + y1 + + a4*hill_promote(A, k4, n4)')
net.add_rate_rule('B', 'a2*hill_repress(C, k2, n2) - b2*B + y2')
net.add_rate_rule('C', 'a3*hill_repress(A, k3, n3) - b3*C + y3')

network = net
int_time = (0, 200)
