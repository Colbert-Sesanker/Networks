% ensemble is a struct of posterior samples and plot options

function plotTraces(ensemble)

param_idxs  = ensemble.traces.params;   

for param_idx = param_idxs
    paramName = ensemble.paramMap{param_idx};    
    figure();
    plot(ensemble.samples(:, param_idx)); 
    xlabel('posterior sample number');
    ylabel(paramName);
    title(['Trace plot for ' paramName]);
end 

end 