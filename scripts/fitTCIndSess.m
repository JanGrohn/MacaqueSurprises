function [betas,names] = fitTCIndSess(DT,region,window)
% average the timecourses for 1, 2, and 3 drops

DT = getEpochs('d1',region, DT, window);

betasRew0 = nan(size(DT.rewEpochs,2),max(DT.monkey),15,6);
names = strings(max(DT.monkey),15,6);


betasRew1 = betasRew0;
betasRew2 = betasRew0;
betasRew3 = betasRew0;


for Asub = 1:max(DT.monkey)
    AsubData = DT(DT.monkey==Asub,:);
    sessions = unique(AsubData.session)';
    for sess = sessions
        SessData = AsubData(AsubData.session==sess,:);
        for epoch = 1:size(SessData.rewEpochs,2)

            display(['fitting monkey ',num2str(Asub),' session ',num2str(sess),'/',num2str(max(AsubData.session)),' epoch ',num2str(epoch)])
            
            SessData.TrialNoNew = nan(height(SessData),1);
            SessData.TrialNoNew(SessData.reward~=0) = 1:length(SessData.TrialNoNew(SessData.reward~=0));
            
            SessData.epochRew = SessData.rewEpochs(:,epoch);
            

            betasRew1(epoch,Asub,sess,:)=nanmean(SessData.epochRew(SessData.reward==1));
            betasRew2(epoch,Asub,sess,:)=nanmean(SessData.epochRew(SessData.reward==2));
            betasRew3(epoch,Asub,sess,:)=nanmean(SessData.epochRew(SessData.reward==3));


        end
    end
end

betas = {};

betas.betasRew1 = betasRew1;
betas.betasRew2 = betasRew2;
betas.betasRew3 = betasRew3;



