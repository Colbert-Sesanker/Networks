%%% Regression of Complex Gene Regulation with using LassoLARS in Octave
%%% 
%%%   Examples:
%%% 1. [a b] = Regress_LassoO('synthmod_01_01.csv', 13, [5,1,6], 'k_ranges', [2,1,5], 'spline_step', .5);
%%% 2. [a b] = Regress_LassoO('synthmod_02_01.csv', 10, [2,1,3], 'k_ranges', [100,10,110]);
%%% 3. [a b] = Regress_LassoO('synthmod_03_01.csv', 10, [2,1,3], 'k_ranges', [100,10,110]);
%%% 4. [a b] = Regress_LassoO('synthmod_04_01_01.csv', 10, [2,1,4], 'k_ranges', [10,20,50]);
%%% 5. [a b] = Regress_LassoO('synthmod_04_02_01.csv', 10, [2,1,4], 'k_ranges', [10,20,50]);
%%% 6. [a b] = Regress_LassoO('synthmod_05_01_01.csv', 50, [2,1,4], 'k_ranges', [300,100,400], 'spline_step', 3);
%%% 7. [a b] = Regress_LassoO('synthmod_05_02_01.csv', 50, [2,1,4], 'k_ranges', [300,100,400],'spline_step', 3);
%%% 
%%%
%%% Required Files: LassoLARS.m, process_options.m (snaged to give Matlab varargin functionality), matrix_sort.m, find_k.m
%%%
%%% Colbert Sesanker Oct 16/2012

function [simple_explanations, complex_explanations] = Regress_LassoO(data, norm_bound, n_ranges, varargin) 

[k_ranges, spline_step] = process_options(varargin, 'k_ranges',1, 'spline_step',5); 

data = csvread(data);            
data = data(2:end,:);              % eliminate 1st row, (labels)
time = data(:,1);                  % time vector is first column of data
samples = length(time);             
nodes = data(:,2:end);             % redefine data without time column 
number_nodes = length(nodes(1,:)); % number of genes/species
splined_time = time(1):spline_step:time(end); % reindex time to spline steps
splined_nodes = zeros(length(splined_time), number_nodes);
for i = 1 : number_nodes
 splined_nodes(:,i) = interp1(time,nodes(:,i),splined_time,'spline'); % cubic spline interpolation with [time(end)-time(start)]/spline_step points
end
nodes = splined_nodes;
time = splined_time';

time_diff = diff(time);
explanator_pairs = nchoosek(number_nodes-1,2);
explanation_parameters = 13;
if (explanator_pairs <= 10)
  complex_candidates = explanator_pairs;
else
  complex_candidates = 10;
end
%% Allocate Matricies:
simple_explanations = zeros(number_nodes, number_nodes); 
complex_explanation_bag = zeros(explanation_parameters,1, number_nodes - 2, number_nodes - 2, number_nodes);
complex_explanations = zeros(explanation_parameters, complex_candidates,number_nodes);
complex_confidence = zeros(number_nodes, number_nodes, number_nodes);

%% Search for Explanations:
for i = 1 : number_nodes      % targets    
  target = nodes(:,i);        % the node being examined for explanation  
  target_diff = diff(target); % vector of differences between x0 and x1 of length(n/2)
  target_derivative = target_diff ./ time_diff; % time derivative of target
  best_parameters = zeros(7,1);
  candidate_parameters = zeros(5,1);      
  one_vector = ones(length(target) - 1, 1);
  prev = 0; 
  for j = 1 : number_nodes     % 1st sources. Note: source1 is never = source2
    if (j ~= i)          
      source1 = nodes(:,j);
      for r = 1 : number_nodes % 2nd sources 
        if (r ~= i && r > j)  % sources pairs are unordered
          best_parameters_error = inf;    % arbitrary large value   
          source2 = nodes(:,r);          
          if length(k_ranges) ~= 1,
            k1s = k_ranges(1,1): k_ranges(1,2) : k_ranges(1,3);
            k2s = k1s;           
          else
            k = find_k(target, source1, source2, spline_step);  
            k1s = k(1,:);
            k2s = k(2,:);
          endif
          for k1 = k1s
            for n1 = n_ranges(1): n_ranges(2): n_ranges(3)                      
                for k2 = k2s
                  for n2 = n_ranges(1): n_ranges(2): n_ranges(3)       
                      hill_function1 = (source1.^n1) ./ (k1^n1 + (source1.^n1));    % activation hill function for source1         
                      hill_function2 = (source2.^n2) ./ (k2^n2 + (source2.^n2));                    
                      hill_function12 = hill_function1 .* hill_function2;  % cross term hill function                                                   
                      design_matrix = [hill_function1(1:end-1), hill_function2(1:end-1), hill_function12(1:end-1)...
                      target(1:end-1), one_vector];
                      % Note: residual_L1 = residual*[1 + L1_norm(candidate_parameters)]   
                      % where residual = L2_norm(design_matrix*candidate_parameters - y)                    
                      [candidate_parameters, residual_L1] = LassoLARS(design_matrix, target_derivative, norm_bound);                        
                      if (residual_L1 < best_parameters_error)
                        best_parameters = [candidate_parameters; [j,r,n1,n2,k1,k2]']; % Best parameters over k's and n's
                        best_parameters_error = residual_L1;                       
                      endif                    
                  endfor 
                endfor
            endfor
          endfor
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
          endif          
         endif
        w = norm(best_parameters(1:3),1); % L1 norm of alphas
        complex_confidence(j,r,i) = best_parameters_error;   % Best parameters over the nCk(number_nodes-1,2) possible pairwise explanations  
        residual = best_parameters_error / (1 + w);
        complex_explanation_bag(:,1,j,r,i) = [best_parameters; [residual, w]']; 
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
end
endfunction
