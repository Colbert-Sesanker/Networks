function plot3dParamSlices(sampledParamEnsemble, figStart)                  
    
figure(figStart);
    
for p = 1:size(sampledParamEnsemble, 1)
    hold on;
    paramTriple = sampledParamEnsemble(p, :);
    x = paramTriple(1);
    y = paramTriple(2);
    z = paramTriple(3);
    scatter3(x, y, z); 
    drawnow;  
      
end

% hXlabel = xlabel('k_1');
% hYlabel = ylabel('n_1');
% hZlabel = zlabel('a_1'); 
hTitle  = title('3d correlation plot params k_1 vs. n_1 vs. a_1');        


set([hXlabel, hYlabel, hZlabel, hTitle], 'FontName', 'AvantGarde');

set([hXlabel, hYlabel, hZlabel], 'FontSize',   12, ...
                                 'FontWeight', 'bold');
set(hTitle,                      'FontSize',   13, ...
                                 'FontWeight','bold');
       
end
