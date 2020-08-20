function plotROIzStat(region,DT,regName,axRange)


regressors = {'sRPE','VS','RRE'};

figure('Renderer', 'painters', 'Position', [50 50 300 400])

C = brewermap(6,'Spectral');

hold on

dist = 0:10:length(regressors)*10;

p = nan(3,4);
chi2 = nan(3,4);

tbl = table;
[VS idx] = unique(DT.(region)(:,2),'stable');
idx = idx(~isnan(VS));
VS = VS(~isnan(VS));
tbl.VS = VS;
tbl.RRE = unique(DT.(region)(~isnan(DT.(region)(:,1)),3),'stable');
tbl.SRPE = unique(DT.(region)(~isnan(DT.(region)(:,1)),1),'stable');
tbl.monkey = DT.monkey(idx);

SRPE = fitlme(tbl,'SRPE~1+(1|monkey)');
RRE = fitlme(tbl,'RRE~1+(1|monkey)');
VS = fitlme(tbl,'VS~1+(1|monkey)');


for reg = 1:3
    if reg == 1
        m = SRPE.Coefficients.Estimate;
        se = SRPE.Coefficients.SE;
    elseif reg == 2
        m = VS.Coefficients.Estimate;
        se = VS.Coefficients.SE;
    elseif reg == 3
        m = RRE.Coefficients.Estimate;
        se = RRE.Coefficients.SE;
    end
    bar(mean([0.75+dist(reg) 6.25+dist(reg)]),m,6,'facecolor', [0.5 0.5 0.5],'edgecolor','none','linewidth',2)
    
    for Asub = 1:6
        
        
     
        for sess = unique(DT.session(DT.monkey==Asub&~isnan(DT.(region)(:,1))))'
                if sum(DT.ChangingDown(DT.monkey==Asub&DT.session==sess)) > 0
                    marker = 'o';
                elseif sum(DT.ChangingUp(DT.monkey==Asub&DT.session==sess)) > 0
                    marker = 'o';
                elseif sum(DT.Equiprobable(DT.monkey==Asub&DT.session==sess)) > 0                
                    marker = 'o';
                elseif sum(DT.Stable(DT.monkey==Asub&DT.session==sess)) > 0
                    marker = 'o';
                else
                    Asub
                    sess
                    error('something is wrong')
                end
                s = scatter(1.75+Asub/2+dist(reg)+randn*0.05,unique(DT.(region)(DT.monkey==Asub&DT.session==sess,reg)),15,C(Asub,:),marker,'filled','markeredgecolor','k');
                s.MarkerFaceAlpha = 0.5;
                s.MarkerEdgeAlpha = 0.3;
            end
    end
        errorbar(mean([0.75+dist(reg) 6.25+dist(reg)]),m,se,'.k','capsize',0,'linewidth',2)
    
end

plot([0 length(regressors)*10],[0 0],'k')

ylabel('z-statistic (a.u.)','fontsize',16)
xticks([3.5:10:length(regressors)*10+3.5])
xticklabels(regressors)
set(gca,'TickLabelInterpreter','none')
xtickangle(45)
ax = gca;
ax.FontSize = 16; 
axis(axRange)
title(regName)
