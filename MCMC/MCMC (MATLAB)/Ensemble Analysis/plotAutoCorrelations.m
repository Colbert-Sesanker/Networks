function plotAutoCorrelations(ensemble)
param_idxs        = ensemble.autoCorr.params;
   
for i = 1: length(param_idxs)    
    paramNames{i} =  ensemble.paramMap{param_idxs(i)}; 
end

maxLags           = ensemble.autoCorr.maxLags;

for  i = 1: length(param_idxs)         
    param_samples = ensemble.samples(:, i);    
    auto_corr     = autocorr(param_samples, maxLags);
    figure();
    plot(auto_corr); 
    xlabel('lag');    
    title(['Auto-Correlation: ' paramNames{i}]);
end 

end 
