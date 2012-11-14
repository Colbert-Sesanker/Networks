%%% Regression of Complex Gene Regulation with LassoLARS in Matlab
%%% 
%%%   Examples:
%%% 1. Regress_LassoM('synthmod_01_01.csv', [2,1,5], [5,1,6]);
%%% 2. Regress_LassoM('synthmod_02_01.csv', [100,10,110], [1,1,3]);
%%% 3. Regress_LassoM('synthmod_03_01.csv', [100,10,110], [1,1,3]);


function [simple_explanations, complex_explanations, best_parameters] = Regress_LassoM(data, k_ranges, n_ranges) 
data = csvread(data, 1);          % read in csv file without 1st row, (labels)
time = data(:,1);                 % time vector is first column of data
nodes = data(:,2:end);            % redefine data without time column 

if size(k_ranges, 1) == 1,
    k_ranges = [k_ranges; k_ranges];
end

if size(n_ranges, 1) == 1,
    n_ranges = [n_ranges; n_ranges];
end

number_nodes = length(nodes(1,:));% number of genes/species
time_diff = diff(time);
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
          for k1 = k_ranges(1,1): k_ranges(1,2) : k_ranges(1,3)
            for n1 = n_ranges(1,1): n_ranges(1,2): n_ranges(1,3)                      
                for k2 = k_ranges(2,1): k_ranges(2,2) : k_ranges(2,3)
                  for n2 = n_ranges(2,1): n_ranges(2,2): n_ranges(2,3)       
                      hill_function1 = (source1.^n1) ./ (k1^n1 + (source1.^n1));    % activation hill function for source1         
                      hill_function2 = (source2.^n2) ./ (k2^n2 + (source2.^n2));                    
                      hill_function12 = hill_function1 .* hill_function2;  % cross term hill function                                                   
                      design_matrix = [hill_function1(1:end-1), hill_function2(1:end-1), hill_function12(1:end-1)...
                      target(1:end-1), one_vector];                        
                      [candidate_parameters, residual] = LassoLARS(design_matrix, target_derivative, 10);                        
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
        beta  = best_parameters(4);         
        new_parameters = (best_parameters_error ~= prev);
         if (beta < 0 && new_parameters)    % reduced dimension simple explanation slice     
          if (abs(alpha1) > abs(alpha2))      
             simple_explanations(i,j) = sign(alpha1)*(1 / best_parameters_error);            
             prev = best_parameters_error;
          else 
             simple_explanations(i,r) = sign(alpha2)*(1 / best_parameters_error);          
             prev = best_parameters_error;
          end          
         end
        complex_confidence(i,j,r) = best_parameters_error;      
        complex_explanations(i,j,r,:) = best_parameters; 
       else
          continue
        end                   
      end        
    else
      continue
    end     
   end
end
end
