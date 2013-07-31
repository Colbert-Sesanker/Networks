function Output  = logNormPDF(Values, Means, Variance)

% Values:     numOfTimePts by 1 Vector
% Means:      numOfTimePts by 1 Vector
% Variance:   numOfTimePts by 1 Matrix

% Make sure values is a column vector
if size(Values, 2) > 1
    Values = Values';
end

% Make sure means is a column vector
if size(Means, 2) > 1
    Means = Means' ;
end

numOfTimePts = length(Values);
assert(length(Means) == length(Values));

Output = sum(  - ones(numOfTimePts, 1) *  (0.5 * log(2*pi*Variance)) -... 
                                                                      ...
                 ((Values - Means).^2) ./ (2*ones(numOfTimePts, 1) * Variance)... 
            );

end
