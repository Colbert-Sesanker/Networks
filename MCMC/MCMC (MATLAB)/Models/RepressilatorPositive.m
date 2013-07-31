function dy = Repressilator_Positive(t, y, p)

%{

 A_initial = .48199;
 B_initial =  5.11385; 
 C_initial =  105.77422; 


 k1 = 2.35804;   k2 = 4.42269;   k3 = 4.80922;   k4 = 5;
 n1 = 5.03005;   n2 = 5.73448;   n3 = 6.05167;   n4 = 7;
 a1 = 5.73702;   a2 = 6.92108;   a3 = 7.46407;   a4 = .7;
 b1 = .3284;     b2 = .4967;     b3 = .4518;
 y1 = .908;      y2 = .8093;     y3 = 1.1444;

%}
numberOfStates     =  3;
numberOfParameters = 18;

%  Species

A   = y(1);
B   = y(2);
C   = y(3);

% Parameters

k1 = p(1);
k2 = p(2);
k3 = p(3);
k4 = p(4);
n1 = p(5);
n2 = p(6);
n3 = p(7);
n4 = p(8);
a1 = p(9);
a2 = p(10);
a3 = p(11);
a4 = p(12);
b1 = p(13);
b2 = p(14);
b3 = p(15);
y1 = p(16);
y2 = p(17);
y3 = p(18);

% Define hill functions
hill_promote = @(x, k, n)  x^n / (k^n + x^n);
hill_repress = @(x, k, n)  k^n / (k^n + x^n);

% Evaluate equations
dy    = zeros(numberOfStates, 1);    % a column vector

dy(1) = a4*hill_promote(A, k4, n4) + a1*hill_repress(B, k1, n1) - b1*A + y1;
dy(2) =                              a2*hill_repress(C, k2, n2) - b2*B + y2;
dy(3) =                              a3*hill_repress(A, k3, n3) - b3*C + y3;


end



