% calculates metric tensor as: 
% G(i,j) = sum_{species}(sens_i(species) * cov_mat * sens_j(species))

function G  = metric_Tensor(...
                             sampledParameters,  Sensitivities_1,... 
                             numSampledParams,   SpeciesObserved,... 
                             CurrentNoise,       prior_second_derivative...                              
                           )
G = zeros(numSampledParams,...
          numSampledParams);    
      
    for SpeciesNum = SpeciesObserved
        for i = 1: numSampledParams
            for j = i: numSampledParams
                G(i, j) = G(i, j) +...                 
                (Sensitivities_1{i}(:, SpeciesNum)' *...
                 Sensitivities_1{j}(:, SpeciesNum)...
                ) ...
                / CurrentNoise(SpeciesNum);                
            end
        end
    end        
    % Expected Fisher information (metric tensor) 
    % is symmetric so set G(i,j) = G(j,i):    
    G    = G + (G - diag(diag(G)))' ;
    
    % Add prior to the metric tensor:    
    % Follows from (log(posterior))_{theta,theta}, second derivative of 
    % log(posterior). Note: this is a prior on a matrix of all sampled parameters        
    G    = G -  prior_second_derivative(numSampledParams,...
                                        sampledParameters);
    % bound singular values near zero
    % cutOff is the maximum condition number allowed
    % 'singVal_cutOff' is a lower bound on the singular values that
    % enforces an upper bound of 'cutOff' on the condition number    
    identity               = eye(numSampledParams);      
    cutOff                 = 1e8; 
    [U, singularValues, V] = svd(G);
    spectralRadius         = max(diag(singularValues));
    singVal_cutOff         = spectralRadius / cutOff;
    bounded_singularValues = max(singularValues, ...
                                 identity  * singVal_cutOff);
    
    G           = U * bounded_singularValues * V';   
    
    % In addition to bounding singular values
    % add diagonal dust to improve rank     
    inv_G       = identity / (G + identity*1e-6);
    
  
end % function   



