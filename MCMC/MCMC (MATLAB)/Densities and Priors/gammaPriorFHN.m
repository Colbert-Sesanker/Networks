function priorProb = gammaPriorFHN(paramNum, param)

%%%%%%%%%
% log Gamma %
%%%%%%%%%
shape = 1;
scale = 3;

if (param < 0)
   priorProb = - Inf;
else
% Calculate probability of value from the prior
   priorProb = log(gampdf(param, shape, scale));

end
      
        
end

