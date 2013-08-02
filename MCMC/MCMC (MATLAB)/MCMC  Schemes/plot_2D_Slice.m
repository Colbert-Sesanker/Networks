function plot_2D_Slice( paramPairs, ...
                        figureNums...
                       )

for p = 1: length(paramPairs)           
        figure(figureNums(p));
        hold on;
        x = paramPairs{p}(1);
        y = paramPairs{p}(2);
        plot(x, y, '-o'); 
        drawnow;
        xlabel(num2str(p));
        %ylabel(num2str(y));
        title(['2d slice pair ' num2str(p)...
               '1 vs 2']);
end
 
end
