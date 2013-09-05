% plots quantile trajectories for specified quantiles 

function plotTrajectoryQuantiles(ensemble)

states             =  ensemble.trajectoryQuantiles.statesToPlot;
numQuantiles       =  length(ensemble.trajQuantiles);

for i = 1: length(states)    
    stateNames{i}  =  ensemble.stateMap{states(i)}; 
end

for i = 1: length(states)
    quantile_trajs = trajectoryQuantiles(ensemble, states(i));

    for j = 1: numQuantiles    
        quantile_trajs{j}
        hold on;
        figure(50 + i);
        plot(quantile_trajs{j}); 
        xlabel('time'); 
        ylabel(stateNames{i});
        title(['95% confidence intervals on trajectories: '...
                stateNames{i}]);
        
    end % for

end % for

end % function
