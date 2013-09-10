function priorProb = gammaPriorDeriv2(paramNum, param)

% Return the partial derivative of the log(prior) w.r.t. the parameter

% log Gamma Prior
shape = 1e7;
scale = 5.737e-7;    

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


