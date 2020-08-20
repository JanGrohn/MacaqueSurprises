function plotFitResult(GLME,GLMEmonkey)

figure('Renderer', 'painters', 'Position', [10 10 400 500])
hold on;

C = brewermap(6,'Spectral');
plot([0 length(GLME.Coefficients)],[0 0],'k')

for Asub = 1:6
    s = scatter([1:length(GLME.Coefficients)-1]-(Asub-3.5)/50,GLMEmonkey{Asub}.Coefficients(2:end,2),30,'filled','markerfacecolor',C(Asub,:),'markeredgecolor','k');
    s.MarkerFaceAlpha = 0.3;
    s.MarkerEdgeAlpha = 0.3;
end
errorbar(GLME.Coefficients(2:end,2),GLME.Coefficients(2:end,3),'ok','linewidth',1,'markersize',5,'capsize',0,'markerfacecolor','k');
for coeff = 2:length(GLME.Coefficients)
    p = double(GLME.Coefficients(coeff,6));
    y = double(GLME.Coefficients(coeff,2)) + double(GLME.Coefficients(coeff,3)) + 0.01;
%     if p < 0.001
%         text(coeff-1,y,'***','horizontalalignment','center','FontSize',16)
%     elseif p < 0.01
%         text(coeff-1,y,'**','horizontalalignment','center','FontSize',16)
%     elseif p < 0.05
%         text(coeff-1,y,'*','horizontalalignment','center','FontSize',16)
%     
%     end
end
            
xticks(1:length(GLME.Coefficients)-1);
xticklabels(GLME.CoefficientNames(2:end));
xtickangle(45)
ylabel('regression weight (a.u.)','FontSize',16)

ax = gca;
ax.FontSize = 16; 

