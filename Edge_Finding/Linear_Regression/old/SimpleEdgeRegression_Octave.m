%%% First read csv file into matrix and extract nodes. 
%%% Chop off labels & time
%%% Edge regression Octave version, search over ks and ns then fit parameters
%%% Colbert Sesanker, 09/2012


data = csvread("synthmod_03_01.csv");      % read in csv file
data = data(2:end,:);              % chop off first row of data
time = data(:,1);                  % time vector is first column of data
time_pts = length(time);           % number of time points
nodes = data(:,2:end);             % redefine data without time column 
                                   % normalize and spline data?
number_nodes = length(nodes(1,:)); % number of genes/species 

%%% 'correlations' is matrix of matricies. The ith matrix offers explanations of node i
%%% in terms of correlations to functions (max, min ...) evaluated pairwise on
%%% the other nodes

time_diff = diff(time);
k_range = 100;
n_range = 3;
explanations = zeros(number_nodes, number_nodes - 1);
confidence = zeros(3,1);
% best_parameters_error = 0
error = 0;
% compute correlations with the other individual nodes

for i = 1 : number_nodes      % targets
  target = nodes(:,i);        % the node being examined for explanation  
  target_diff = diff(target); % vector of differences between x0 and x1 of length(n/2)
  target_derivative = target_diff ./ time_diff; % Computes differences as first derivatives
  best_parameters = zeros(3,1);
  previous_error = 1000;
  for j = 1 : number_nodes    % sources
      if (j ~= i)          
        source = nodes(:,j);        
        for k = 1:  k_range
            for n = 1:  n_range
              hill_function = (source.^n) ./ (k^n + (source.^n)); 
              one_vector = ones(length(target) - 1, 1);
              design_matrix = [hill_function(1:end - 1), target(1:end - 1), one_vector]; 
              [candidate_parameters,_, error] = LinearRegression(design_matrix, target_derivative);
              if (error < previous_error)
                 best_parameters = candidate_parameters;
                 best_parameters_error = error;
              endif
              previous_error = error;
            endfor 
        endfor
      else
        continue
      endif
      alpha = best_parameters(1);
      beta  = best_parameters(2);
      gamma = best_parameters(3);
      confidence(i,j) = best_parameters_error;
      if (alpha > 0 && beta < 0)        
        explanations(i,j) = 1;
      elseif (alpha < 0 && beta < 0)  
        explanations(i,j) = -1;     
      endif     
   endfor
endfor
               
