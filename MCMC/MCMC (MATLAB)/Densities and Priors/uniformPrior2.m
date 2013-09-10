function priorProb = uniformPrior2(paramNum, param)

%%%%%%%%%%%%%%%
% log Uniform %
%%%%%%%%%%%%%%%

a = 5.72;
b = 5.754;

% Calculate probability of value from the prior
if paramNum == 9             
   priorProb = log(unifpdf(param, a, b));
else    
   priorProb = 0;
end      
        
end

