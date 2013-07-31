function dy = FitzHughNagumo(t, y, p)
 
% Initial conditions for FHN model are V = -1, R = 1

% Set up species
V  = y(1);
R  = y(2);


% Set up parameters
a = p(1);
b = p(2);
c = p(3);

% Evaluate equations
dy    = zeros(2,1);    % a column vector

dy(1) = c*(V - (V^3) / 3 + R);
dy(2) = - (V - a + b * R) / c;

end
