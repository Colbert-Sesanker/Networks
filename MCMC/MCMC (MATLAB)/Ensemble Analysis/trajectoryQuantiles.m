% returns quantile trajectories for specified quantiles 

function quantile_trajs = trajectoryQuantiles(ensemble, state)

quantiles         = ensemble.trajQuantiles;
numQuantiles      = length(quantiles);
 
% Trajectories are rows, timePoints are columns
% want a submatix here. should be a cleaner way to do this
trajectories(:,:)  = ensemble.trajectories(state, :, :); 

% At each time point, sort the trajectories, this 
% will yeild a set of trajectory qualtiles, none of which are actual
% trajectories
sortedTrajs       = sort(trajectories);
% Number of trajectories 
numTrajs          = size(sortedTrajs, 1);
quantile_trajs    = cell(1, numQuantiles);

for i  = 1: numQuantiles    
    
    q     = quantiles(i);
    index = q * numTrajs;
    
    if isinteger(index)
        quantile_trajs{i} = sortedTrajs(:, index);
    else
        floor_idx = floor(index);
        ceil_idx  = ceil(index);     
        
        quantile_trajs{i} = sortedTrajs(floor_idx, :) + ((index - floor_idx) / (ceil_idx - floor_idx)) *...
                                                                                                        ...
                                                        (sortedTrajs(ceil_idx, :) - sortedTrajs(floor_idx, :));
    end % if   
    
end % for

end 
