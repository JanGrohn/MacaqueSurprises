function plotFitResultCoeff(GLMEDrift,GLMESimple,GLMEEqui,GLMEDriftMonkey,GLMESimpleMonkey,GLMEEquiMonkey,coeff)

figure('Renderer', 'painters', 'Position', [10 10 400 500])
hold on;

C = brewermap(6,'Spectral');
plot([0 4],[0 0],'k')

for Asub = 1:6
    s = scatter(1-(Asub-3.5)/50,GLMEDriftMonkey{Asub}.Coefficients(coeff,2),30,'filled','markerfacecolor',C(Asub,:),'markeredgecolor','k');
%     text(1-(Asub-3.5)/50,GLMEDriftMonkey{Asub}.Coefficients(coeff,2),num2str(Asub));
        s.MarkerFaceAlpha = 0.3;
    s.MarkerEdgeAlpha = 0.3;
    s = scatter(2-(Asub-3.5)/50,GLMESimpleMonkey{Asub}.Coefficients(coeff,2),30,'filled','markerfacecolor',C(Asub,:),'markeredgecolor','k');
% text(2-(Asub-3.5)/50,GLMESimpleMonkey{Asub}.Coefficients(coeff,2),num2str(Asub));
        s.MarkerFaceAlpha = 0.3;
    s.MarkerEdgeAlpha = 0.3;
    s = scatter(3-(Asub-3.5)/50,GLMEEquiMonkey{Asub}.Coefficients(coeff,2),30,'filled','markerfacecolor',C(Asub,:),'markeredgecolor','k');
% text(3-(Asub-3.5)/50,GLMEEquiMonkey{Asub}.Coefficients(coeff,2),num2str(Asub));
        s.MarkerFaceAlpha = 0.3;
    s.MarkerEdgeAlpha = 0.3;
end
errorbar(1,GLMEDrift.Coefficients(coeff,2),GLMEDrift.Coefficients(coeff,3),'ok','linewidth',1,'markersize',5,'capsize',0,'markerfacecolor','k');
errorbar(2,GLMESimple.Coefficients(coeff,2),GLMESimple.Coefficients(coeff,3),'ok','linewidth',1,'markersize',5,'capsize',0,'markerfacecolor','k');
errorbar(3,GLMEEqui.Coefficients(coeff,2),GLMEEqui.Coefficients(coeff,3),'ok','linewidth',1,'markersize',5,'capsize',0,'markerfacecolor','k');

           
xticks(1:3);
xticklabels({'RRE (Volatile)','RRE (Stable)','RRE (Equiprobable)'});
xtickangle(45)
ylabel('regression weight (a.u.)','FontSize',16)

ax = gca;
ax.FontSize = 16; 

