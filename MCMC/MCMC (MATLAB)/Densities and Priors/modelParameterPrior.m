function PP = ModelParameterPrior(ParaNum, param)

%%%%%%%%%
% log Gamma %
%%%%%%%%%

if (param < 0)
        PP = - Inf;
else
% Calculate probability of value from the prior
PP = log(gampdf(param, 1, 3));
end
      
        
end

