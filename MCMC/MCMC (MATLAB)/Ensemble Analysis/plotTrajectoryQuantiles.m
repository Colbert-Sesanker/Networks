% plots quantile trajectories for specified quantiles 

function quantile_trajs = plotTrajectoryQuantiles(ensemble, state)

stateName      =  ensemble.stateMap{state}; 
quantile_trajs = trajectoryQuantiles(ensemble, state);

for i = 1: length(quantile_trajs)
    
    hold on;
    figure(50);
    plot(quantile_trajs{i}); 
    xlabel('time'); 
    ylabel(stateName);
    title(['95% confidence intervals on trajectories: ' stateName]);
    
end % for

end 
