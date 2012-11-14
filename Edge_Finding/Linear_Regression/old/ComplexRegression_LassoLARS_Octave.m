%%% First read csv file into matrix and extract nodes. 
%%% Chop off labels & time
%%% Complex Regulation LassoLARS Octave


data = csvread("synthmod_01_01.csv");      % read in csv file
data = data(2:end,:);             % chop off first row of data
time = data(:,1);                 % time vector is first column of data
time_pts = length(time);          % number of time points
nodes = data(:,2:end);            % redefine data without time column 
                                 
number_nodes = length(nodes(1,:));% number of genes/species 

time_diff = diff(time);
k1_range = 5;
n1_range = 6;
k2_range = 5;
n2_range = 6;
simple_explanations = zeros(number_nodes, number_nodes); 
complex_explanations = zeros(number_nodes, number_nodes - 1, number_nodes - 1,5);
complex_confidence = zeros(number_nodes, number_nodes, number_nodes);

for i = 1 : number_nodes      % targets
  target = nodes(:,i);        % the node being examined for explanation  
  target_diff = diff(target); % vector of differences between x0 and x1 of length(n/2)
  target_derivative = target_diff ./ time_diff; % time derivative of target
  best_parameters = zeros(5,1);
  candidate_parameters = zeros(5,1);
  best_parameters_error = 0;
  previous_error = 1000;      % arbitrary large value
  one_vector = ones(length(target) - 1, 1);
  prev = 0; 
  for j = 1 : number_nodes    % 1st sources. Note: source1 is never = source2
    if (j ~= i)          
      source1 = nodes(:,j);
      for r = 1 : number_nodes % 2nd sources 
        if (r ~= i && r > j)  % sources pairs are unorderedbest_parameters_error
          source2 = nodes(:,r);
          for k1 = 1: k1_range
            for n1 = 3:  n1_range                      
                for k2 = 1: k2_range
                  for n2 = 3: n2_range     
                      hill_function1 = (source1.^n1) ./ (k1^n1 + (source1.^n1));    % activation hill function for source1         
                      hill_function2 = (source2.^n2) ./ (k2^n2 + (source2.^n2));                    
                      hill_function12 = hill_function1 .* hill_function2;  % cross term hill function                                                   
                      design_matrix = [hill_function1(1:end-1), hill_function2(1:end-1), hill_function12(1:end-1)...
                      target(1:end-1), one_vector];                        
                      [candidate_parameters, residual] = LassoShooting2(design_matrix, target_derivative, 1);% 10);                        
                      if (residual < previous_error)
                        best_parameters = candidate_parameters;
                        best_parameters_error = residual;                       
                      endif
                      previous_error = residual;
                  endfor 
                endfor
            endfor
          endfor
        alpha1 = best_parameters(1);
        alpha2 = best_parameters(2);
        alpha12 = best_parameters(3);     
        beta  = best_parameters(4);
        gamma = best_parameters(5);               
        same = sign(alpha1) == sign(alpha2);
        new_parameters = (best_parameters_error ~= prev);
         if (beta < 0 && new_parameters)    % reduced dimension simple explanation slice     
          if (abs(alpha1) > abs(alpha2))      
             simple_explanations(i,j) = sign(alpha1)*(1 / best_parameters_error);            
             prev = best_parameters_error;
          elseif 
             simple_explanations(i,r) = sign(alpha2)*(1 / best_parameters_error);          
             prev = best_parameters_error;
          endif          
        endif
        complex_confidence(i,j,r) = best_parameters_error;      
        complex_explanations(i,j,r,:) = best_parameters; 
       else
          continue
        endif                   
      endfor        
    else
      continue
    endif      
   endfor
endfor               
