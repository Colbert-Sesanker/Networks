function dy = network_fiveDotTwo(t, y, p)



numberOfStates     =  5;

%  Species

A   = y(1);
B   = y(2);
C   = y(3);
D   = y(4);
E   = y(5);


% Parameters

k1 = p(1);
k2 = p(2);
k3 = p(3);
k4 = p(4);
k5 = p(5);
k6 = p(6);
k7 = p(7);

n1 = p(8);
n2 = p(9);
n3 = p(10);
n4 = p(11);
n5 = p(12);
n6 = p(13);
n7 = p(14);


a1 = p(15);
a2 = p(16);
a3 = p(17);
a4 = p(18);
a5 = p(19);

b1 = p(20);
b2 = p(21);
b3 = p(22);
b4 = p(23);
b5 = p(24);

y1 = p(25);
y2 = p(26);
y3 = p(27);
y4 = p(28);
y5 = p(29);


% Evaluate equations
dy    = zeros(numberOfStates, 1);    % a column vector

dy(1) = a1*(k1^n1/(k1^n1 + D^n1))                          - b1*A + y1;
dy(2) = a2*(A^n2 / (k2^n2 + A^n2))                         - b2*B + y2;

dy(3) = a3*(A^n3 / (k3^n3 + A^n3))*(B^n4 / (k4^n4 + B^n4)) - b3*C + y3;
dy(4) = a4*(B^n5 / (k5^n5 + B^n5))*(C^n6 / (k6^n6 + C^n6)) - b4*D + y4;

dy(5) = a5*(A^n7 / (k7^n7 + A^n7))                         - b5*E + y5; 

end % function



