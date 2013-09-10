function priorProb = gammaPrior(paramNum, param)

%%%%%%%%%%%%%
% log Gamma %
%%%%%%%%%%%%%

if (param < 0)
   priorProb = - Inf;
else
% Calculate probability of value from the prior
    if paramNum == 9
       shape = 1;
       scale = 6;       
       priorProb = log(gampdf(param, shape, scale));
    else    
       priorProb = 0;
    end      
        
end

