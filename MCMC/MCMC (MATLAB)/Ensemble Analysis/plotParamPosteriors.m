function plotParamPosteriors(ensemble)
posteriorParams   = ensemble.samples;   
param_idxs        = ensemble.posteriorParamsToPlot;             
bins              = ensemble.posteriorParamPlotBins;
axisRanges        = ensemble.posteriorParamAxisRanges;
                    
for p     = 1:length(param_idxs)
    param_idx     = param_idxs(p);
    paramName     = ensemble.paramMap{param_idx}; 
    param_samples = posteriorParams(:, p);      
    figure(100 + p);
    histfit(param_samples, bins);
    axis(axisRanges{p});
    hXlabel = xlabel('param value');    
    hTitle  = title(['Posterior of Param ' paramName]);
    set([hXlabel, hTitle], 'FontName', 'AvantGarde');
    set(hXlabel, 'FontSize',  12,...
                 'FontWeight', 'bold');
    set(hTitle, 'FontSize',   13, ...
                'FontWeight', 'bold');
end 

end 
