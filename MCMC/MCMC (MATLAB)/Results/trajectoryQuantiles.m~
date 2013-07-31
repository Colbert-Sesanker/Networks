% In 'trajectories', Trajectories are columns, timePoints are rows
function quantile_trajs = trajectoryQuantiles(trajectories,...
                                              quantiles)

numQuantiles      = length(quantiles);
sortedTrajs       = sort(trajectories);
timePoints        = rows(sortedTrajs);
quantile_trajs    = cell(numQuantiles);

for i  = 1: numQuantiles    
    q     = quantiles(i);
    index = q * timePoints;
    
    if isinteger(index)
        quantile_trajs{i} = sortedTrajs(:, index);
    else
        floor_idx = floor(index);
        ceil_idx  = ceil(index);
        q_floor   = floor_idx / timePoints;
        q_ceil    = ceil_idx  / timePoints;
        
        quantile_trajs{i} = sortedTrajs(floor_idx, :) + ((q - q_floor) / (q_ceil - q_floor)) *...
                                                                                              ...
                                                        (sortedTrajs(ceil_idx, :) - sortedTrajs(floor_idx, :));
   end % if
    
    
end % for

end % function
