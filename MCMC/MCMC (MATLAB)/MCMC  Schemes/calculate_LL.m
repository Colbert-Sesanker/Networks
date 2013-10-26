% Calculates log-likelihood of estimated data
% LL is for log-likelihood

function LL = calculate_LL(speciesEstimates, Y,... 
                           currentNoise, speciesObserved) 
                       
  for n = speciesObserved         
          
      LL(n) = logNormPDF(speciesEstimates(n, :), ...
                         Y(n, :), currentNoise(n));
  end
  
  % This sum is in [-inf, 0]
  LL = sum(LL);
   
end

