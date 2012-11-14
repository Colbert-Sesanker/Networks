%%% First read csv file into matrix and extract nodes. 
%%% Chop off labels & time

data = csvread('synthmod_02_01.csv', 1); % read csv excluding first row, the labels 
time = data(:,1);                 % time vector is first column of data
time_pts = length(time);          % number of time points
nodes = data(:,2:end);             % redefine data without time column 
                                  % normalize and spline data?
number_nodes = length(nodes(1,:)); % take length of arbitrary column

time_diff = diff(time);
k1_range = 10;
n1_range = 10;
k2_range = 10;
n2_range = 10;
simple_explanations = zeros(number_nodes, number_nodes); 
simple_confidence = zeros(number_nodes, number_nodes);
complex_explanations = zeros(number_nodes, number_nodes - 1, number_nodes - 1);
complex_confidence = zeros(number_nodes, number_nodes, number_nodes);


for i = 1 : number_nodes      % targets
  target = nodes(:,i);        % the node being examined for explanation  
  target_diff = diff(target); % vector of differences between x0 and x1 of length(n/2)
  target_derivative = target_diff ./ time_diff; % time derivative of target
  best_parameters = zeros(5,1);
  candidate_parameters = zeros(5,1);
  best_parameters_error = 0;
  previous_error = 1000; % arbitrary large value
  residual = 0;          % error
  for j = 1 : number_nodes    % 1st sources
    if (j ~= i)          
      source1 = nodes(:,j);
      for r = 1 : number_nodes % 2nd sources 
        if (r ~= i && r > j)  % sources pairs are unordered
          source2 = nodes(:,r);
          for k1 = 1:  k1_range
            for n1 = 1:  n1_range
              hill_function1 = (source1.^n1) ./ (k1^n1 + (source1.^n1));    % activation hill function for source1
                for k2 = 1: k2_range
                  for n2 =1: n2_range             
                      hill_function2 = (source2.^n2) ./ (k2^n2 + (source2.^n2));   
                      hill_function12 = hill_function1 .* hill_function2;  % cross term hill function
                      one_vector = ones(length(target) - 1, 1);
                      design_matrix = [hill_function1(1:end-1), hill_function2(1:end-1), hill_function12(1:end-1)...
                      target(1:end-1), one_vector]; % create matrix from coulmn vectors TRANSPOSE HILL FUNC in matlab
                      [candidate_parameters,na, residual] = regress(target_derivative, design_matrix);
                      residual = norm(residual, 2);
                      if (residual < previous_error)
                        best_parameters = candidate_parameters;
                        best_parameters_error = residual;
                      end
                      previous_error = residual;
                  end 
                end
            end
          end
        alpha1 = best_parameters(1);
        alpha2 = best_parameters(2);
        alpha12 = best_parameters(3);     
        beta  = best_parameters(4);
        gamma = best_parameters(5);               
        same = sign(alpha1) == sign(alpha2);
        if (beta < 0)
          if (abs(alpha1) > abs(alpha2))              
             simple_explanations(i,j) = sign(alpha1);
             simple_confidence(i,j) = best_parameters_error;
          else
             simple_explanations(i,r) = sign(alpha2);
             simple_confidence(i,r) = best_parameters_error;
          end          
        end
        complex_confidence(i,j,r) = best_parameters_error;        
        if (abs(alpha12) > abs(alpha1 + alpha2))
          if (alpha12 > 0)
            complex_explanations(i,j,r) = 2;            
          else 
            complex_explanations(i,j,r) = -2;
            
          end                  
        elseif (same && abs(alpha12) < abs(alpha1 + alpha2))
          if (alpha1 > 0)
            complex_explanations(i,j,r) = 1;
          else
            complex_explanations(i,j,r) = -1;
          end
        end          
        else
          continue
        end                   
      end
    else
      continue
    end
      
   end
end
               
