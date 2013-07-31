function prior_prob = ModelParameterLogPriorDerivative(ParaNum, param)

% Return the partial derivative of the log(prior) w.r.t. the parameter

% log Gamma Prior
a = 1;
b = 3;

if (param < 0)
     error(['Parameter cannot be negative:'...
            'redoing current iteration'])                   
            
else
    prior_prob = (a - 1) / param - (1 / b);
end


end


