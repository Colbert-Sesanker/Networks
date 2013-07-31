function plotParameterQuantiles(ensemble)
param_idxs        = ensemble.paramQuantiles.params; 

for param_idx     = param_idxs;
    paramName     = ensemble.paramMap{param_idx}; 
    param_samples = ensemble.samples(:, param_idx);    
    
    figure();
    hist(param_samples);
    xlabel('parameter value'); 
    ylabel('frequency');
    title(['Posterior Density: ' paramName]);
end 

end 
