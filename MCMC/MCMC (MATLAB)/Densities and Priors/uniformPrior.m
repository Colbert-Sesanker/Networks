function priorProb = uniformPrior(paramNum,              param,...
                                  paramIdxsForThisPrior, interval)
                              
% The indexes here are the absolute indexes for the parameter

%%%%%%%%%%%%%%%
% log Uniform %
%%%%%%%%%%%%%%%

a = interval(1);
b = interval(2);

% Calculate probability of value from the prior
if any(paramNum == paramIdxsForThisPrior)              
   priorProb = log(unifpdf(param, a, b));
else    
   priorProb = 0;
end      
        
end

