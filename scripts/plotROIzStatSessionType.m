function plotROIzStatSessionType(region,DT,regName,axRange)


regressors = {'sRPE','VS','RRE'};

figure('Renderer', 'painters', 'Position', [50 50 600 600])

C = brewermap(6,'Spectral');

hold on

plot([0 length(regressors)*15],[0 0],'k')
dist = 0:15:length(regressors)*15;

tbl = table;
[VS idx] = unique(DT.(region)(:,2),'stable');
idx = idx(~isnan(VS));
VS = VS(~isnan(VS));
tbl.VS = VS;
tbl.RRE = unique(DT.(region)(~isnan(DT.(region)(:,1)),3),'stable');
tbl.SRPE = unique(DT.(region)(~isnan(DT.(region)(:,1)),1),'stable');
tbl.monkey = DT.monkey(idx);
tbl.Changing = DT.Changing(idx);
tbl.Stable = DT.Stable(idx);
tbl.Equi = DT.Equiprobable(idx);


SRPElea = fitlme(tbl(tbl.Changing==1,:),'SRPE~1+(1|monkey)');
SRPEsta = fitlme(tbl(tbl.Stable==1,:),'SRPE~1+(1|monkey)');
SRPEequi = fitlme(tbl(tbl.Equi==1,:),'SRPE~1+(1|monkey)');
RRElea = fitlme(tbl(tbl.Changing==1,:),'RRE~1+(1|monkey)');
RREsta = fitlme(tbl(tbl.Stable==1,:),'RRE~1+(1|monkey)');
RREequi = fitlme(tbl(tbl.Equi==1,:),'RRE~1+(1|monkey)');
VSlea = fitlme(tbl(tbl.Changing==1,:),'VS~1+(1|monkey)');
VSsta = fitlme(tbl(tbl.Stable==1,:),'VS~1+(1|monkey)');
VSequi = fitlme(tbl(tbl.Equi==1,:),'VS~1+(1|monkey)');



for reg = 1:3


    if reg == 1
        m_lea = SRPElea.Coefficients.Estimate;
        se_lea = SRPElea.Coefficients.SE;
        m_sta = SRPEsta.Coefficients.Estimate;
        se_sta = SRPEsta.Coefficients.SE;
        m_equi = SRPEequi.Coefficients.Estimate;
        se_equi = SRPEequi.Coefficients.SE;
    elseif reg == 3
        m_lea = RRElea.Coefficients.Estimate;
        se_lea = RRElea.Coefficients.SE;
        m_sta = RREsta.Coefficients.Estimate;
        se_sta = RREsta.Coefficients.SE;
        m_equi = RREequi.Coefficients.Estimate;
        se_equi = RREequi.Coefficients.SE;
    elseif reg == 2
        m_lea = VSlea.Coefficients.Estimate;
        se_lea = VSlea.Coefficients.SE;
        m_sta = VSsta.Coefficients.Estimate;
        se_sta = VSsta.Coefficients.SE;
        m_equi = VSequi.Coefficients.Estimate;
        se_equi = VSequi.Coefficients.SE;
    end
    bar(mean([0.25+dist(reg) 3.25+dist(reg)]),m_lea,3,'facecolor', [0.5 0.5 0.5],'edgecolor','none','linewidth',2)
    bar(mean([4.25+dist(reg) 7.25+dist(reg)]),m_sta,3,'facecolor', [0.5 0.5 0.5],'edgecolor','none','linewidth',2)
    bar(mean([8.25+dist(reg) 11.25+dist(reg)]),m_equi,3,'facecolor', [0.5 0.5 0.5],'edgecolor','none','linewidth',2)
    
    for Asub = 1:6
        for sess = unique(DT.session(DT.monkey==Asub&~isnan(DT.(region)(:,1))))'
                if sum(DT.ChangingDown(DT.monkey==Asub&DT.session==sess)) > 0
                    marker = 'o';
                    x = Asub/4+dist(reg)+randn*0.05;
                elseif sum(DT.ChangingUp(DT.monkey==Asub&DT.session==sess)) > 0
                    marker = 'o';
                    x = Asub/4+dist(reg)+randn*0.05;
                elseif sum(DT.Equiprobable(DT.monkey==Asub&DT.session==sess)) > 0                
                    marker = 'o';
                    x = 8+Asub/4+dist(reg)+randn*0.05;
                elseif sum(DT.Stable(DT.monkey==Asub&DT.session==sess)) > 0
                    marker = 'o';
                    x = 4+Asub/4+dist(reg)+randn*0.05;
                else
                    Asub
                    sess
                    error('something is wrong')
                end
                s = scatter(.75+x,unique(DT.(region)(DT.monkey==Asub&DT.session==sess,reg)),15,C(Asub,:),marker,'filled','markeredgecolor','k');
                s.MarkerFaceAlpha = 0.5;
                s.MarkerEdgeAlpha = 0.3;
            end
        end
    means = nan(6,1);
    ses = nan(6,1);
    meansVol = nan(6,1);
    meansStable = nan(6,1);
    meansEqui = nan(6,1);
    sesVol = nan(6,1);
    sesStable = nan(6,1);
    sesEqui = nan(6,1);
    for Asub = 1:6             
        means(Asub) = mean(unique(DT.(region)(DT.monkey==Asub&~isnan(DT.(region)(:,1)),reg)));
        ses(Asub) = std(unique(DT.(region)(DT.monkey==Asub&~isnan(DT.(region)(:,1)),reg)))/sqrt(length(unique(DT.(region)(DT.monkey==Asub&~isnan(DT.(region)(:,1)),reg))));   
        
        meansVol(Asub) = mean(unique(DT.(region)(DT.Changing==1&DT.monkey==Asub&~isnan(DT.(region)(:,1)),reg)));
        meansStable(Asub) = mean(unique(DT.(region)(DT.Stable==1&DT.monkey==Asub&~isnan(DT.(region)(:,1)),reg)));
        meansEqui(Asub) = mean(unique(DT.(region)(DT.Equiprobable==1&DT.monkey==Asub&~isnan(DT.(region)(:,1)),reg)));
        
        sesVol(Asub) = std(unique(DT.(region)(DT.Changing==1&DT.monkey==Asub&~isnan(DT.(region)(:,1)),reg)))/sqrt(length(unique(DT.(region)(DT.Changing==1&DT.monkey==Asub&~isnan(DT.(region)(:,1)),reg))));
        sesStable(Asub) = std(unique(DT.(region)(DT.Stable==1&DT.monkey==Asub&~isnan(DT.(region)(:,1)),reg)))/sqrt(length(unique(DT.(region)(DT.Stable==1&DT.monkey==Asub&~isnan(DT.(region)(:,1)),reg))));
        sesEqui(Asub) = std(unique(DT.(region)(DT.Equiprobable==1&DT.monkey==Asub&~isnan(DT.(region)(:,1)),reg)))/sqrt(length(unique(DT.(region)(DT.Equiprobable==1&DT.monkey==Asub&~isnan(DT.(region)(:,1)),reg))));
    end
    



    errorbar(mean([0.25+dist(reg) 3.25+dist(reg)]),m_lea,se_lea,'.k','capsize',0,'linewidth',2)
    errorbar(mean([4.25+dist(reg) 7.25+dist(reg)]),m_sta,se_sta,'.k','capsize',0,'linewidth',2)
    errorbar(mean([8.25+dist(reg) 11.25+dist(reg)]),m_equi,se_equi,'.k','capsize',0,'linewidth',2)

end

regressors = {'sRPE (changing)','sRPE (stable)','sRPE (equiprobable)','VS (changing)','VS (stable)','VS (equiprobable)','RRE (changing)','RRE (stable)','RRE (equiprobable)'};

ylabel('z-statistic (a.u.)','fontsize',16)
xticks([1.75 5.75 9.75 1.75+15 5.75+15 9.75+15 1.75+30 5.75+30 9.75+30])


axis(axRange)
xticklabels(regressors)
set(gca,'TickLabelInterpreter','none')
xtickangle(45)
ax = gca;
ax.FontSize = 16; 
title(regName)