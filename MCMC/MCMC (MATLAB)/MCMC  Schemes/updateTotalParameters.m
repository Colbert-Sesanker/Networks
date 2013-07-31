% Updates total parameters with the new values of sampled parameters
function tp = update_totalParameters( totalParameters,...
                                      newSampledParameters,...
                                      sampled_Param_idxs)
% copy total parameters                                 
tp                     = totalParameters;
% set initial sampling params of tp to new values                                                     
tp(sampled_Param_idxs) = newSampledParameters;
end


