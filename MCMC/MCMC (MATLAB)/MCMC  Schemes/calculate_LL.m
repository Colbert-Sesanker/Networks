% Calculates log-likelihood of estimated data
% LL is for log-likelihood

function LL = calculate_LL(speciesEstimates, Y,... 
                           CurrentNoise,     SpeciesObserved) 
                       
  for n = SpeciesObserved         
          
      LL(n) = logNormPDF(speciesEstimates(n, :), ...
                         Y(n, :), CurrentNoise(n));
  end
  
  LL = sum(LL);
   
end

