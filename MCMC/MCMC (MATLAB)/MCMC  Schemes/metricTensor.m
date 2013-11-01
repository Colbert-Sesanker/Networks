 % calculates metric tensor as: 
% G(i,j) = sum_{species}(sens_i(species) * cov_mat * sens_j(species))

function G  = metricTensor(...
                             sampledParameters,  sensitivities_1,... 
                             numSampledParams,   speciesObserved,... 
                             currentNoise,       priorSecondDerivative...                            
                          )
G = zeros(numSampledParams,...
          numSampledParams);    
      
    for speciesNum = speciesObserved
        for i = 1: numSampledParams
            for j = i: numSampledParams
                G(i, j) = G(i, j) +...                 
                (sensitivities_1{i}(:, speciesNum)' *...
                 sensitivities_1{j}(:, speciesNum)...
                ) ...
                / currentNoise(speciesNum);                
            end
        end
    end        
    % Expected Fisher information (metric tensor) 
    % is symmetric so set G(i,j) = G(j,i):    
    G    = G  + (G - diag(diag(G)))' ;
    
    % Add prior to the metric tensor:    
    % Follows from (log(posterior))_{theta,theta}, second derivative of 
    % log(posterior). Note: this is a prior on a matrix of all sampled parameters        
    G    = G  -  priorSecondDerivative(numSampledParams,...
                                       sampledParameters);
    % bound singular values near zero
    % cutOff is the maximum condition number allowed
    % 'singVal_cutOff' is a lower bound on the singular values that
    % enforces an upper bound of 'cutOff' on the condition number    
    identity               = eye(numSampledParams);      
    cutOff                 = 1e10; 
    [U, singularValues, V] = svd(G);
    spectralRadius         = max(diag(singularValues));
    singVal_cutOff         = spectralRadius / cutOff;
    bounded_singularValues = max(singularValues, ...
                                 identity  * singVal_cutOff);
    
    G           = U * bounded_singularValues * V'; 
    % Add diagonal dust to improve conditioning
    %G           = G + identity*1e-6;
    
    % check if this matrix is dangerous and should raise an error
    % rcond is the reciprocal condition number
    if  hasNaN(G) % ||  rcond(G) < 1e-12 )
        error('has NaN values')
    end    

    
%     % multiply G by beta in [0,1] for tempering
%     beta = 1;   % setting beta == 1 makes G no longer the mertic tensor
%     G = beta*G;
    
    % In addition to bounding singular values
    % add diagonal dust to improve rank     
    %invG       = identity / (G + identity*1e-6);
    
  
end % function   



