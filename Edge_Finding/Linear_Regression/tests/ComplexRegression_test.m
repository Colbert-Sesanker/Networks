%% For testing output of know regression functions

data = csvread("synthmod_03_01.csv");      % read in csv file
data = data(2:end,:);             % chop off first row of data
time = data(:,1);                 % time vector is first column of data
time_pts = length(time);          % number of time points
nodes = data(:,2:end);            % redefine data without time column 
                                  % normalize and spline data?
number_nodes = length(nodes(1,:));% number of genes/species 

time_diff = diff(time);
k1_range = 10;
n1_range = 10;
k2_range = 10;
n2_range = 10;
simple_explanations = zeros(number_nodes, number_nodes); 
simple_confidence = zeros(number_nodes, number_nodes);
complex_explanations = zeros(number_nodes, number_nodes - 1, number_nodes - 1);
complex_confidence = zeros(number_nodes, number_nodes, number_nodes);

%%% 'correlations' is matrix of matricies. The ith matrix offers explanations of node i
%%% in terms of correlations to functions (max, min ...) evaluated pairwise on
%%% the other nodes

time_diff = diff(time);
k_range = 100;
n_range = 5;
previous_error = 1000;
explanations = zeros(number_nodes, number_nodes - 1); 
best_parameters = zeros(3,1);
confidence = best_parameters = zeros(3,1);
error = 0;


target = nodes(:,3);        % the node being examined for explanation  
source = nodes(:,2); 
target_diff = diff(target); % vector of differences between x0 and x1 of length(n/2)
target_derivative = target_diff ./ time_diff;  
k = 100;
n = 3;

hill_function = (source.^n) ./ (k^n + (source.^n));
one_vector = ones(length(target) - 1, 1);
design_matrix = [hill_function(1:end -1), target(1:end-1), one_vector]; % create matrix from coulmn vectors TRANSPOSE HILL FUNC in matlab
[candidate_parameters,_, error] = LinearRegression(design_matrix, target_derivative);
      


dif = zeros(rows(design_matrix),1);  
regressed_y_prime = zeros(length(dif),1);
for i = 1 : rows(design_matrix)
  regressed_y_prime(i) = design_matrix(i,:) * candidate_parameters;  % [10,-.05,.05]';      % candidate_parameters;
  y_prime = target_derivative(i);
  dif(i) = y_prime -regressed_y_prime(i);
end

check = [target_derivative(1:49), regressed_y_prime, dif];

error = norm(dif); 
   

          
