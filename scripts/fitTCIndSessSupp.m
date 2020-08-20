function [betas,names] = fitTCIndSessSupp(DT,region,window)
% average the timecourses for 1 - 10 drops

DT = getEpochs('d2',region, DT, window);

betasRew0 = nan(size(DT.rewEpochs,2),max(DT.monkey),15,9);
names = strings(max(DT.monkey),15,1);

betasRew1 = betasRew0;
betasRew2 = betasRew0;
betasRew3 = betasRew0;
betasRew4 = betasRew0;
betasRew5 = betasRew0;
betasRew6 = betasRew0;
betasRew7 = betasRew0;
betasRew8 = betasRew0;
betasRew9 = betasRew0;
betasRew10 = betasRew0;



for Asub = 1:max(DT.monkey)
    AsubData = DT(DT.monkey==Asub,:);
    sessions = unique(AsubData.session)';
    for sess = sessions
        SessData = AsubData(AsubData.session==sess,:);
        for epoch = 1:size(SessData.rewEpochs,2)

            display(['fitting session ',num2str(sess),'/',num2str(max(AsubData.session)),' epoch ',num2str(epoch)])

            SessData.epochRew = SessData.rewEpochs(:,epoch);
            

            betasRew1(epoch,Asub,sess,:)=nanmean(SessData.epochRew(SessData.outcome==1));
            betasRew2(epoch,Asub,sess,:)=nanmean(SessData.epochRew(SessData.outcome==2));
            betasRew3(epoch,Asub,sess,:)=nanmean(SessData.epochRew(SessData.outcome==3));
            betasRew4(epoch,Asub,sess,:)=nanmean(SessData.epochRew(SessData.outcome==4));
            betasRew5(epoch,Asub,sess,:)=nanmean(SessData.epochRew(SessData.outcome==5));
            betasRew6(epoch,Asub,sess,:)=nanmean(SessData.epochRew(SessData.outcome==6));
            betasRew7(epoch,Asub,sess,:)=nanmean(SessData.epochRew(SessData.outcome==7));
            betasRew8(epoch,Asub,sess,:)=nanmean(SessData.epochRew(SessData.outcome==8));
            betasRew9(epoch,Asub,sess,:)=nanmean(SessData.epochRew(SessData.outcome==9));
            betasRew10(epoch,Asub,sess,:)=nanmean(SessData.epochRew(SessData.outcome==10));

        end
    end
end



names(Asub,sess,:) = "(Intercept)";

betas = {};

betas.betasRew1 = betasRew1;
betas.betasRew2 = betasRew2;
betas.betasRew3 = betasRew3;
betas.betasRew4 = betasRew4;
betas.betasRew5 = betasRew5;
betas.betasRew6 = betasRew6;
betas.betasRew7 = betasRew7;
betas.betasRew8 = betasRew8;
betas.betasRew9 = betasRew9;
betas.betasRew10 = betasRew10;




