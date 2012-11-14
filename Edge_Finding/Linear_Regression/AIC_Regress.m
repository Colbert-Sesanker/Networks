%%%  First read csv file into matrix and extract nodes. 
%%%  Chop off labels & time
%%%  Examples:
%%%  source("AIC_Regress.m")
%%%  [a b] = AIC_Regress('synthmod_04_02_01.csv', [10,20,50], [2,1,4]);
%%%  [a b] = AIC_Regress('synthmod_05_01_01.csv', [300,100,400], [2,1,4]);



function [simple_explanations, complex_explanations] = AIC_Regress(data, k_ranges, n_ranges) 
data = csvread(data);             % read in csv file without 1st row, (labels)
data = data(2:end,:);   
time = data(:,1);                 % time vector is first column of data
samples = length(time);
nodes = data(:,2:end);            % redefine data without time column 

if size(k_ranges, 1) == 1,
    k_ranges = [k_ranges; k_ranges];
end

if size(n_ranges, 1) == 1,
    n_ranges = [n_ranges; n_ranges];
end

number_nodes = length(nodes(1,:)); % number of genes/species
time_diff = diff(time);
explanator_pairs = nchoosek(number_nodes-1,2);
explanation_parameters = 9;
if (explanator_pairs <= 10)
  complex_candidates = explanator_pairs;
else
  complex_candidates = 10;
end
simple_explanations = zeros(number_nodes, number_nodes); 
complex_explanation_bag = zeros(explanation_parameters,1, number_nodes - 2, number_nodes - 2, number_nodes);
complex_explanations = zeros(explanation_parameters, complex_candidates,number_nodes);
complex_confidence = zeros(number_nodes, number_nodes, number_nodes);

for i = 1 : number_nodes      % targets    
  target = nodes(:,i);        % the node being examined for explanation  
  target_diff = diff(target); % vector of differences between x0 and x1 of length(n/2)
  target_derivative = target_diff ./ time_diff; % time derivative of target
  best_parameters = zeros(7,1);
  candidate_parameters = zeros(5,1);
  best_parameters_error = 0;
  previous_error = 1000;      % arbitrary large value
  one_vector = ones(length(target) - 1, 1);
  prev = 0; 
  for j = 1 : number_nodes     % 1st sources. Note: source1 is never = source2
    if (j ~= i)          
      source1 = nodes(:,j);
      for r = 1 : number_nodes % 2nd sources 
        if (r ~= i && r > j)  % sources pairs are unorderedbest_parameters_error
          source2 = nodes(:,r);
          for k1 = k_ranges(1,1): k_ranges(1,2) : k_ranges(1,3)
            for n1 = n_ranges(1,1): n_ranges(1,2): n_ranges(1,3)                    
              hill_function1 = (source1.^n1) ./ (k1^n1 + (source1.^n1));    % activation hill function for source1
              design_matrix = [hill_function1(1:end-1), target(1:end-1), one_vector]; 
              [parameters,_, residual(1)] = LinearRegression(design_matrix, target_derivative); 
              candidate_parameters(1,1) = parameters(1);
              candidate_parameters(4,1) = parameters(2);
              candidate_parameters(5,1) = parameters(3);
                 for k2 = k_ranges(2,1): k_ranges(2,2) : k_ranges(2,3)
                  for n2 = n_ranges(2,1): n_ranges(2,2): n_ranges(2,3)            
                      hill_function2 = (source2.^n2) ./ (k2^n2 + (source2.^n2));
                      design_matrix = [hill_function2(1:end-1), target(1:end-1), one_vector];
                      [parameters,_, residual(2)] = LinearRegression(design_matrix, target_derivative);  
                      candidate_parameters(2,2) = parameters(1);
                      candidate_parameters(4,2) = parameters(2);
                      candidate_parameters(5,2) = parameters(3);
                   
                      design_matrix = [hill_function1(1:end-1), hill_function2(1:end-1), target(1:end-1), one_vector];
                      [parameters,_, residual(3)] = LinearRegression(design_matrix, target_derivative); 
                      candidate_parameters(1,3) = parameters(1);
                      candidate_parameters(2,3) = parameters(2);
                      candidate_parameters(4,3) = parameters(3);
                      candidate_parameters(5,3) = parameters(4);                        
                      
                      hill_function12 = hill_function1 .* hill_function2;  % cross term hill function
                 
                      design_matrix = [hill_function12(1:end-1),target(1:end-1), one_vector];
                      [parameters,_, residual(4)] = LinearRegression(design_matrix, target_derivative); 
                      candidate_parameters(3,4) = parameters(1);
                      candidate_parameters(4,4) = parameters(2);
                      candidate_parameters(5,4) = parameters(3);
                      
                      design_matrix = [hill_function1(1:end-1), hill_function12(1:end-1),target(1:end-1), one_vector];
                      [parameters,_, residual(5)] = LinearRegression(design_matrix, target_derivative); 
                      candidate_parameters(1,5) = parameters(1);
                      candidate_parameters(3,5) = parameters(2);
                      candidate_parameters(4,5) = parameters(3);
                      candidate_parameters(5,5) = parameters(4);     

                      design_matrix = [hill_function2(1:end-1), hill_function12(1:end-1),target(1:end-1), one_vector];
                      [parameters,_, residual(6)] = LinearRegression(design_matrix, target_derivative);  
                      candidate_parameters(2,6) = parameters(1);
                      candidate_parameters(3,6) = parameters(2);
                      candidate_parameters(4,6) = parameters(3);
                      candidate_parameters(5,6) = parameters(4);    
                                           
                      design_matrix = [hill_function1(1:end-1), hill_function2(1:end-1), hill_function12(1:end-1)...
                      target(1:end-1), one_vector];                        
                      [parameters,_, residual(7)] = LinearRegression(design_matrix, target_derivative); 
                      candidate_parameters(1:5,7) = parameters;
                      
                      mr = max(residual); % residual = residual / max(residual);
                      RSS = residual.^2;
                      % Select model based on Akaikie Information Criterion
                      for s= 1: length(residual)
                       bp = candidate_parameters(1:5,s);  
                       k = nnz(bp(1:3));   
                       kn = (k+norm(bp(1:3),1)); 
                       AIC = samples*log(RSS(s)/samples) + 2*kn; 
                       AICc = AIC + 2*kn*(kn+1)/(samples-kn-1);
                       AICc = (nnz(bp)/ length(bp)) + residual(s)/ mr; % Higher penalty, no justification,
                       penalized_residual(s) = AICc;         
                      end
                      [best_residual, best_residual_index] = min(penalized_residual);        
                      if (best_residual < previous_error)
                        best_parameters(:) = [candidate_parameters(:,best_residual_index); [j,r]'];
                        best_parameters_error = penalized_residual(best_residual_index);  
                        best_residual_error = residual(best_residual_index);                     
                      endif
                      previous_error = best_residual;
                  endfor 
                endfor
            endfor
          endfor
        alpha1 = best_parameters(1);
        alpha2 = best_parameters(2);          
        beta  = best_parameters(4);              
        new_parameters = (best_residual_error ~= prev);
         if (beta < 0 && new_parameters)    % reduced dimension simple explanation slice     
          if (abs(alpha1) > abs(alpha2))      
             simple_explanations(i,j) = sign(alpha1)*(1 / best_residual_error);            
             prev = best_residual_error;
          else 
             simple_explanations(i,r) = sign(alpha2)*(1 / best_residual_error);          
             prev = best_residual_error;
          endif          
         endif 
        RSS = best_residual_error^2;       
        k = nnz(best_parameters);       
        AIC = samples*log(RSS/samples) + 2*k;       
        AICc= AIC + 2*k*(k+1)/(samples-k-1);
        complex_confidence(j,r,i) = AICc;   % Best parameters over the nCk(number_nodes-1,2) possible pairwise explanations  
        complex_explanation_bag(:,1,j,r,i) = [best_parameters; [best_residual_error, AICc]']; 
       else
          continue
        endif                  
      endfor        
    else
      continue
    endif     
   endfor
endfor
indexes = zeros(complex_candidates, 2);
for i = 1: number_nodes
  c = complex_confidence(:,:,i);
  c(c==0) = inf; % replace all zeros in matrix with infinity for sorting
  [na, indexes] = matrix_sort(c, complex_candidates)
  for j = 1: complex_candidates
     l = indexes(j,1)
     r = indexes(j,2)
     complex_explanations(:,j,i) = complex_explanation_bag(:,1,l,r,i);
  endfor
% complex_explanation_bag(i,[indexes(:,1)],[indexes(:,2)],:)
end
endfunction
