%%% First read csv file into matrix and extract nodes. 
%%% Chop off labels & time

data = csvread("synthmod_01_01.csv");      % read in csv file
data = data(2:end,:);             % chop off first row of data
time = data(:,1);                 % time vector is first column of data
time_pts = length(time);          % number of time points
nodes = data(:,2:end);            % redefine data without time column 
                                  % normalize and spline data?
number_nodes = length(nodes(1,:));% number of genes/species 

time_diff = diff(time);
k1_range = 5;
n1_range = 6;
k2_range = 5;
n2_range = 6;
simple_explanations = zeros(number_nodes, number_nodes); 
complex_explanations = zeros(number_nodes, number_nodes - 1, number_nodes - 1);
complex_confidence = zeros(number_nodes, number_nodes, number_nodes);
%best = zeros(5, number_nodes);

for i = 1 : number_nodes      % targets
  target = nodes(:,i);        % the node being examined for explanation  
  target_diff = diff(target); % vector of differences between x0 and x1 of length(n/2)
  target_derivative = target_diff ./ time_diff; % time derivative of target
  best_parameters = zeros(5,1);
  candidate_parameters = zeros(5,7);
  best_parameters_error = 0;
  previous_error = 1000; % arbitrary large value
  residual = zeros(1,7);          % error
  one_vector = ones(length(target) - 1, 1);
  j_old = 1;
  r_old = 1;
  for j = 1 : number_nodes    % 1st sources. Note: source1 is never = source2
    if (j ~= i)          
      source1 = nodes(:,j);
      for r = 1 : number_nodes % 2nd sources 
        if (r ~= i && r > j)  % sources pairs are unordered
          source2 = nodes(:,r);
          for k1 = 2: k1_range
            for n1 = 5:  n1_range

              hill_function1 = (source1.^n1) ./ (k1^n1 + (source1.^n1));    % activation hill function for source1

              design_matrix = [hill_function1(1:end-1), target(1:end-1), one_vector]; 
              [parameters,_, residual(1)] = LinearRegression(design_matrix, target_derivative); 
              candidate_parameters(1,1) = parameters(1);
              candidate_parameters(4,1) = parameters(2);
              candidate_parameters(5,1) = parameters(3);

                for k2 = 2: k2_range
                  for n2 = 5: n2_range    
       
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
                      candidate_parameters(:,7) = parameters;
                      
                      mr = max(residual); % residual = residual / max(residual);
                      sample_size = time_pts - 1;
                      RSS = residual.^2;
                      for s= 1: length(residual)
                       bp = candidate_parameters(:,s);  
                       k = nnz(bp);       
                       AIC = sample_size*log(RSS(s)/sample_size) + 2*k;       
                       AICc= AIC + 2*k*(k+1)/(sample_size-k-1);
                       penalized_residual(s) = AICc; %2*(nnz(bp)/ length(bp)) + residual(s)/ mr; 
                      end
                      [best_residual, best_residual_index] = min(penalized_residual);                  
                    
                      %penalized_residual = ((nnz(bp)^20) / length(bp))*best_residual;                 

                      if (best_residual < previous_error)
                        best_parameters = candidate_parameters(:,best_residual_index);
                        best_parameters_error = residual(best_residual_index);                       
                      endif
                      previous_error = best_residual;   
                   
#{     
                      mr = max(residual); % residual = residual / max(residual);
                      [best_residual, best_residual_index] = min(residual);  
                      bp = candidate_parameters(:,best_residual_index);
                      penalized_residual = 5*(nnz(bp)/ length(bp)) + best_residual/mr; 
                      %penalized_residual = ((nnz(bp)^20) / length(bp))*best_residual;                 

                      if (penalized_residual < previous_error)
                        best_parameters = candidate_parameters(:,best_residual_index);
                        best_parameters_error = best_residual;                       
                      endif
                      previous_error = penalized_residual;
#}
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
         if (beta < 0)
          if (abs(alpha1) > abs(alpha2) && simple_explanations(i,j_old) != sign(alpha1)*(1 / best_parameters_error))              
             simple_explanations(i,j) = sign(alpha1)*(1 / best_parameters_error);            
             j_old = j;
          elseif (simple_explanations(i,r_old) != sign(alpha2)*(1 / best_parameters_error))
             simple_explanations(i,r) = sign(alpha2)*(1 / best_parameters_error);          
             r_old = r;
          endif          
        endif
        complex_confidence(i,j,r) = best_parameters_error;        
        if (abs(alpha12) > abs(alpha1 + alpha2))
          if (alpha12 > 0)
            complex_explanations(i,j,r) = 2;            
          else
            complex_explanations(i,j,r) = -2;            
          endif              
        elseif (same && abs(alpha12) < abs(alpha1 + alpha2))
          if (alpha1 > 0)
            complex_explanations(i,j,r) = 1;
          else
            complex_explanations(i,j,r) = -1;
          endif
        endif         
        else
          continue
        endif                   
      endfor        
    else
      continue
    endif      
   endfor
endfor               
