function dy = FitzHughNagumoSens2(t,y,p)

% Initial conditions for FHN model are V = -1, R = 1

% Set up species
V  = y(1);
R  = y(2);
Va = y(3);
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
a = p(1);
b = p(2);
c = p(3);

% Evaluate equations
dy    = zeros(2,1);    % a column vector

dy(1) = c*(V-(V^3)/3+R);
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
