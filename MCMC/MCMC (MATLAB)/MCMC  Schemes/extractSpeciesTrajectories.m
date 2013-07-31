% Extract species time series from trajectories 
function speciesEstimates = extract_Species_Trajectories(trajectories,...
                                                         NumOfSpecies)                                      
   % Note the transpose. Species are rows, times are cols  
  speciesEstimates  = trajectories(:, 1: NumOfSpecies)';   
end
