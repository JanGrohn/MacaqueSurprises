function plotTC(betas,reg,names,plottitle,firstxlabel,DT,window,dat)


t = -1:2.28/10:window;

expNames = fieldnames(betas);
if dat == "d1"
    C = brewermap(length(expNames),'Dark2');
elseif dat == "d2"
    C = brewermap(10,'Blues');
end


mins = nan(length(expNames),1);
maxs = nan(length(expNames),1);



for exp = 1:length(expNames)
    plot([-20 -19],[0 0],'color',C(exp,:))
end

plot([-10 window],[0 0],'k')
plot([0 0],[-10 window],'k')
for exp = 1:length(expNames)

    tc = nan(size(betas.(expNames{exp}),2),length(t));
    for ID = 1:size(betas.(expNames{exp}),2)
        tc(ID,:) = nanmean(betas.(expNames{exp})(:,ID,:,squeeze(names{1}(1,1,:)==reg)),3);
      
    end
    mins(exp) = min(mean(tc));
    maxs(exp) = max(mean(tc));
    shadedErrorBar(t,mean(tc),std(tc)/sqrt(size(betas.(expNames{exp}),2)),{'color',C(exp,:)},1)

    axis([-1 window min(mins)-0.05 max(maxs)+0.05])
    xlabel('time (s)')
    ylabel('mean signal change (a.u.)')
    if exp==1
        title(plottitle)
    end
    xticks([0,5,10])
    xticklabels([firstxlabel,"5","10"])
end

ax = gca;
ax.FontSize = 16; 



