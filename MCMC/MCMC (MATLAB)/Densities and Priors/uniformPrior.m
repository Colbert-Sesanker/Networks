function priorProb = uniformPrior(paramNum, param)

%%%%%%%%%%%%%%%
% log Uniform %
%%%%%%%%%%%%%%%

a = 3;
b = 9;

% Calculate probability of value from the prior
if paramNum == 9             
   priorProb = log(unifpdf(param, a, b));
else    
   priorProb = 0;
end      
        
end

