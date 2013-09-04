function plotAutoCorrelations(ensemble)
param_idxs        = ensemble.autoCorr.params;
   
for i = 1: length(param_idxs)    
    paramNames{i} =  ensemble.paramMap{param_idxs(i)}; 
end

maxLags           = ensemble.autoCorr.maxLags;

for  i = 1: length(param_idxs)   

    param_idx     = param_idxs(i);
    paramName     = ensemble.paramMap{param_idx}; 
    param_samples = ensemble.samples(:, param_idx);    
    auto_corr     = autocorr(param_samples, maxLags);
    figure();
    plot(auto_corr); 
    xlabel('lag');    
    title(['Auto-Correlation: ' paramNames{}]);
end 

end 
