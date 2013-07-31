% Get sensitivities using finite differences
% Returns a cell array of matrices. The columns of each matrix 
% holds the trajectory of (X_p)' for some species X and parameter p 
% species estimates is provided so it doesn't need to be recomputed

function [Sens1, Sens2] = ...
                          ...
         getSensitivities_FD(  equations,          numStates,... 
                               numSampledParams,   totalParameters,...
                               sampledParam_idxs,  speciesEstimatesBest,...
                               timePoints,         initialValues,...
                               epsilon,            zeroMetricTensorDerivatives...     
                            )
                             
    
    sensitivities_1  = cell(1, numSampledParams);
    sensitivities_2  = cell(numSampledParams, numSampledParams);
    
    % this function handle solves ODEs for incremented parameters
    integrate_states = @(incrementedParams)... 
                       integrateEquations( equations,     timePoints,...
                                           initialValues, incrementedParams...
                                          );
                                          
    % initialize incremented parameters to parameters  
    f                =  speciesEstimatesBest' ;   % note transpose
    params           =  totalParameters(sampledParam_idxs);
    % increment and decrement all sampled parameters, simultaneously by (param*epsilon)   
     
    for i = 1 :  numSampledParams
    % Get first order sensitivities of all species with respect to parameter j
    % first central difference: (f(x + epsilon, y) - f(x - epsilon, y)) / 2*epsilon 
    % Note: all f occurrences are matrices
         
         % reset incremented params
         params_inc            =  totalParameters; % reset incremented params
         params_dec            =  totalParameters;
        
         param_idx             =  sampledParam_idxs(i);
         param                 =  params(i);
          
         params_inc(param_idx) =  param  + epsilon*param;                                        
         params_dec(param_idx) =  param  - epsilon*param; 
         
              
         [ ~ , f_inc]          =  integrate_states(params_inc);         
         [ ~ , f_dec]          =  integrate_states(params_dec);               
        
         central_1st           =  (f_inc    -    f_dec) ./ 2*epsilon;
         central_2nd           =  (f_inc - 2*f + f_dec) ./ epsilon^2;                               
        
         sensitivities_1{i}    =  central_1st;
         sensitivities_2{i, i} =  central_2nd;
         
         f_Inc{i} = f_inc; 
         f_Dec{i} = f_dec;
    end
         
    % Get second order sensitivities of all species with respect to parameters i and j
    if zeroMetricTensorDerivatives % zero metric tensor derivatives is false
       
        sensitivities_2 = NaN;
     
     else   
        
        for i =          1: numSampledParams 
            for j =  i + 1: numSampledParams
                params_inc_2              =  totalParameters; % reset incremented params
                params_dec_2              =  totalParameters;
             
                params_idxs               =  sampledParam_idxs([ i j ]);
                param_s                   =  params([ i j ]);       
              
                params_inc_2(params_idxs) =  param_s  + epsilon*param_s;
                params_dec_2(params_idxs) =  param_s  - epsilon*param_s; 
                
                % f_inc_2 = f(x + eps, y + eps), and f_dec_2 = f(x - eps, y - eps)
                [ ~ , f_inc_2]            = integrate_states(params_inc_2);
                [ ~ , f_dec_2]            = integrate_states(params_dec_2);   
             
                central_mixed_2nd         = ( f_inc_2 ...
                                                     - f_Inc{i} - f_Inc{j} + ...
                                              2*f... 
                                                     - f_Dec{i} - f_Dec{j} + ...
                                              f_dec_2...
                                            )... 
                                                ./ 2*epsilon^2;
                       
                sensitivities_2{i, j}     = central_mixed_2nd;                                               
                sensitivities_2{j, i}     = sensitivities_2{i, j};            
                
            end % for
        
        end % for
    end % if
                
                Sens1                     = sensitivities_1;
                Sens2                     = sensitivities_2;
end % function


