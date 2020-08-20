function [betas,names] = fitTCIndSessAlpha(DT,region,window)
% run a regression against the timecourse with sRPE, RRE, VS, and TrialNo
% as regressors. sRPE is constructed using different alphas

DT = getEpochs('d1',region(1), DT, window);
DT2 = getEpochs('d1',region(2), DT, window);
DT.rewEpochs = (DT.rewEpochs+DT2.rewEpochs)/2;

betasRew0 = nan(size(DT.rewEpochs,2),max(DT.monkey),15,5);
names = strings(max(DT.monkey),15,5);

betasRew = betasRew0;

tstatRew = betasRew;

modelstring = 'epochRew~sRPEtemp+RREcurrent+VS+TrialNo';

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
           
            SessData.sRPEtemp = zscore(SessData.sRPEtemp);

            SessData.RREcurrent = zscore(SessData.RREcurrent);
            SessData.VS = zscore(SessData.VS);
            SessData.TrialNo = zscore(SessData.TrialNo);

            fit = fitlm(SessData,modelstring);
            betasRew(epoch,Asub,sess,:) = table2array(fit.Coefficients(:,1));
            tstatRew(epoch,Asub,sess,:) = table2array(fit.Coefficients(:,3));
            names(Asub,sess,:) = string(fit.CoefficientNames);


        end
    end
end

betas = {};


betas.betasRew = betasRew;





