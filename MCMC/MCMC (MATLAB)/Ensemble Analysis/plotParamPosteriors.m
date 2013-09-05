function plotParamPosteriors(ensemble)
posteriorParams   = ensemble.samples;   
param_idxs        = ensemble.posteriorParamsToPlot;             
bins              = ensemble.posteriorParamPlotBins;

for p     = 1:length(param_idxs)
    param_idx     = param_idxs(p);
    paramName     = ensemble.paramMap{param_idx}; 
    param_samples = posteriorParams(:, p);      
    figure(100 + p);
    hist(param_samples, bins);
    xlabel('param value');    
    title(['Posterior of Param ' paramName]);
end 

end 
