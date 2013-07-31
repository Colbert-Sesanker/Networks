function plotTraces(params,    ...
                    figureNums,...
                    iterationNum)

for p = 1: length(params)           
        figure(figureNums(p));
        hold on;
        plot(iterationNum, params(p), '-o'); 
        xlabel('sample number');
        ylabel(num2str(p));
        title(['Trace plot for param ' num2str(p)]);
end
 
end
