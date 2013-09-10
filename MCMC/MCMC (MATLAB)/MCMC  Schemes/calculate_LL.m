% Calculates log-likelihood of estimated data
% LL is for log-likelihood

function LL = calculate_LL(speciesEstimates, Y,... 
                           currentNoise, speciesObserved) 
                       
  for n = speciesObserved         
          
      LL(n) = logNormPDF(speciesEstimates(n, :), ...
                         Y(n, :), currentNoise(n));
  end
  
  LL = sum(LL);
   
end

