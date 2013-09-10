function priorProb = gammaPriorDeriv(paramNum, param)

% Return the partial derivative of the log(prior) w.r.t. the parameter

% log Gamma Prior
shape = 1;
scale = 6;

if (param < 0)
     error(['Parameter cannot be negative:'...
            'redoing current iteration'])                   
            
else
    if paramNum == 9
        priorProb = (shape - 1) / param - (1 / scale);
    else 
        priorProb = 0;
    end
    
end


end


