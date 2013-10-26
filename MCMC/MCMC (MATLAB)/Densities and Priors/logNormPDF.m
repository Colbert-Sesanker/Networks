function output  = logNormPDF(values, means, variance)

% Values:     numOfTimePts by 1 Vector
% Means:      numOfTimePts by 1 Vector
% Variance:   numOfTimePts by 1 Matrix

% Make sure 'values' is a column vector
if size(values, 2) > 1
    values = values';
end

% Make sure means is a column vector
if size(means, 2) > 1
    means = means' ;
end

numOfTimePts = length(values);
assert(length(means) == length(values));

% Assumes isotropic variance
% log of product of 1-D Gaussians for every point of solution trajectory
% for a single state variable. 
% This sum and each of its components are in [-inf, 0]
output = sum(  - ones(numOfTimePts, 1) *  (0.5 * log(2*pi*variance)) -... 
                                                                      ...
                 ((values - means).^2) ./ (2*ones(numOfTimePts, 1) * variance)... 
            );

end
