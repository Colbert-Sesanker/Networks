function priorProb = uninformedSupportPrior(paramNum,               param,...
                                            paramIdxsForThisPrior,  interval)
                              
% The indexes here are the absolute indexes for the parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% probability of 1 on support   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = interval(1);
b = interval(2);

% Calculate probability of value from the prior
if any(paramNum == paramIdxsForThisPrior) 
   % Check if param is in support region
   if param >= a && param <= b 
       priorProb = 0;     % log(1)
   else
       priorProb = -inf;  % log(0)
   end
else    
   priorProb = 0;
end      
        
end

