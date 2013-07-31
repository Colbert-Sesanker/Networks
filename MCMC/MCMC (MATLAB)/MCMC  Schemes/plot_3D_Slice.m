function plot_3D_Slice( paramTriples, ...
                        figureNums...
                       )

for p = 1: length(paramTriples)           
        figure(figureNums(p));
        hold on;
        x = paramTriples{p}(1);
        y = paramTriples{p}(2);
        z = paramTriples{p}(3);
        scatter3(x, y, z); 
        drawnow;
        xlabel(num2str(p));
        %ylabel(num2str(y));
        title(['3d slice params ' num2str(x)...
                num2str(y) num2str(z)]);
end
 
end
