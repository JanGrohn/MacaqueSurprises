%% Data and path
addpath_all;

load('data.mat') % behavioural and neural data (behaviour contains 5 additional session for which we don't have fMRI data)
load('suppdata.mat') % data for Figure S1


%% Fig 2A
% run behavioural models for all animals combined
GLME1 = fitglme(DT,'Lapse~VS+sRE+RRE+Position+TrialNo+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');
GLME1nosRE = fitglme(DT,'Lapse~VS+RRE+Position+TrialNo+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');
GLME1noVS = fitglme(DT,'Lapse~sRE+RRE+Position+TrialNo+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');
GLME1noRRE = fitglme(DT,'Lapse~sRE+VS+Position+TrialNo+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');
GLME1noTrialNo = fitglme(DT,'Lapse~VS+sRE+RRE+Position+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');

% and for each monkey separately 
GLME1monkey = {};
for Asub = 1:6
    GLME1monkey{Asub} = fitglme(DT(DT.monkey==Asub,:),'Lapse~VS+sRE+RRE+Position+TrialNo+(1|monkey:session)','distribution','binomial');
end

% get the statistics by doing likelihood ratio tests
compare(GLME1nosRE,GLME1)
compare(GLME1noVS,GLME1)
compare(GLME1noRRE,GLME1)
compare(GLME1noTrialNo,GLME1)

% plot results
plotFitResult(GLME1,GLME1monkey)


%% Fig 2B
% run behavioural models for all animals combined, this time with separately regressors for the reward history
GLME2 = fitglme(DT,'Lapse~VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');
GLME2not1 = fitglme(DT,'Lapse~VS+t2+t3+t4+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');
GLME2not2 = fitglme(DT,'Lapse~VS+t1+t3+t4+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');
GLME2not3 = fitglme(DT,'Lapse~VS+t1+t2+t4+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');
GLME2not4 = fitglme(DT,'Lapse~VS+t1+t2+t3+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');
GLME2not5 = fitglme(DT,'Lapse~VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');
GLME2noVS = fitglme(DT,'Lapse~t1+t2+t3+t4+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');
GLME2noRRE = fitglme(DT,'Lapse~VS+t1+t2+t3+t4+t5+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');
GLME2noTrialNo = fitglme(DT,'Lapse~VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');

% for each monkey separately 
GLME2monkey = {};
slopes = nan(6,1);
for Asub = 1:6
    GLME2monkey{Asub} = fitglme(DT(DT.monkey==Asub,:),'Lapse~VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo+(1|monkey:session)','distribution','binomial');
    % also fit a line through the reward history of individual monkeys
    p = polyfit(GLME2monkey{Asub}.Coefficients.Estimate(4:8)',1:5,1);
    slopes(Asub) = p(2);
end

% test whether the slope of the line is different from 0
[~, p, ~, s] = ttest(slopes)

% plot the results
plotFitResult(GLME2,GLME2monkey)


%% Fig 3C
% plot the z-statistics extracted from the whole-brain in an ROI placed at
% the peak activity for the sRPE regressor
plotROIzStat('PE_Cluster2',DTn,'right striatum',[0 27 -1.4 2.6])
plotROIzStat('PE_Cluster1',DTn,'left striatum',[0 27 -1.4 2.6])

% get the statistics (for all sessions combined and also for the individual
% session types)
[chi2 p] = getzStatp('PE_Cluster2',DTn)
[chi2 p] = getzStatp('PE_Cluster1',DTn)


%% Fig 3D
% plot the z-statistics extracted from the whole-brain in the a priori
% VTA/SN ROI.
plotROIzStat('SN_r',DTn,'right VTA/SN',[0 27 -2.6 2.2])
plotROIzStat('SN_l',DTn,'left VTA/SN',[0 27 -2.6 2.2])

% get the statistics (for all sessions combined and also for the individual
% session types)
[chi2 p] = getzStatp('SN_r',DTn)
[chi2 p] = getzStatp('SN_l',DTn)
% get the statistics when combining the left and right VTA/SN
[chi2 p] = getzStatp2regions('SN_r','SN_l',DTn)


%% Fig 3E
% get the average BOLD timecourse for 1, 2, and 3 drops in striatum,
% VTA/SN and orofacial cortex
[betas,names,regions] = fitTC(DTn,10,'d1');

% combine the data from the left and right hemispheres
betasAv{3}.betasRew1 = permute([permute(betas{1}.betasRew1,[3 1 2 4]);permute(betas{2}.betasRew1,[3 1 2 4])],[2 3 1 4]);
betasAv{3}.betasRew2 = permute([permute(betas{1}.betasRew2,[3 1 2 4]);permute(betas{2}.betasRew2,[3 1 2 4])],[2 3 1 4]);
betasAv{3}.betasRew3 = permute([permute(betas{1}.betasRew3,[3 1 2 4]);permute(betas{2}.betasRew3,[3 1 2 4])],[2 3 1 4]);

betasAv{1}.betasRew1 = permute([permute(betas{3}.betasRew1,[3 1 2 4]);permute(betas{4}.betasRew1,[3 1 2 4])],[2 3 1 4]);
betasAv{1}.betasRew2 = permute([permute(betas{3}.betasRew2,[3 1 2 4]);permute(betas{4}.betasRew2,[3 1 2 4])],[2 3 1 4]);
betasAv{1}.betasRew3 = permute([permute(betas{3}.betasRew3,[3 1 2 4]);permute(betas{4}.betasRew3,[3 1 2 4])],[2 3 1 4]);

betasAv{2}.betasRew1 = permute([permute(betas{5}.betasRew1,[3 1 2 4]);permute(betas{6}.betasRew1,[3 1 2 4])],[2 3 1 4]);
betasAv{2}.betasRew2 = permute([permute(betas{5}.betasRew2,[3 1 2 4]);permute(betas{6}.betasRew2,[3 1 2 4])],[2 3 1 4]);
betasAv{2}.betasRew3 = permute([permute(betas{5}.betasRew3,[3 1 2 4]);permute(betas{6}.betasRew3,[3 1 2 4])],[2 3 1 4]);

% we plot the Intercept (i.e. the average)
names{1}(1,1,1)='(Intercept)';
reg = '(Intercept)';

plottitle = '';
firstxlabel = "reward delivery";
regions = {'striatum','VTA/SN','orofacial cortex'};

figure;hold on
for ii = 1:length(betasAv)
    subplot(1,3,ii);hold on
    plotTC(betasAv{ii},reg,names,plottitle,firstxlabel,DT,10,"d1")

    title(regions{ii})
end

% get the statistics by computing the area under the curve
stats = integralTC(betasAv,1)
stats = integralTC(betasAv,2)
stats = integralTC(betasAv,3)


%% Fig 4C
% plot the z-statistics extracted from the whole-brain in an ROI placed at
% the peak activity for the RRE regressor
plotROIzStatSessionTypeChangingStable('CurrentTwoDrop_Cluster1',DTn,'striatum',[0 inf -2 1.5])
% get the statistics (for all sessions combined and also for the individual
% session types)
[chi2 p] = getzStatp('CurrentTwoDrop_Cluster1',DTn)


%% Fig 4D
% plot the z-statistics extracted from the whole-brain in an ROI placed at
% the peak activity for the RRE regressor
plotROIzStatSessionTypeChangingStable('CurrentTwoDrop_lOFC',DTn,'plOFC',[0 inf -2 1.5])
% get the statistics (for all sessions combined and also for the individual
% session types)
[chi2 p] = getzStatp('CurrentTwoDrop_lOFC',DTn)


%% Fig 5B
% plot the z-statistics extracted from the whole-brain in an ROI placed at
% the peak activity for the VS regressor
plotROIzStat('SSurprise_Cluster3',DTn,'lPFC',[0 27 -inf 2.3])
% get the statistics (for all sessions combined and also for the individual
% session types)
[chi2 p] = getzStatp('SSurprise_Cluster3',DTn)


%% Fig S1A
% get the average choice (L-R) for each juice amount
av_ch = nan(19,4);
for mag_in=1:19
    mag_dif=mag_in-10;
    for monkey = 1:4
        sub_data=suppDT(suppDT.mag_dif==mag_dif&suppDT.monkey==monkey,:);
        av_ch(mag_in,monkey)=mean(sub_data.choice);
        tri_count(mag_in)=size(sub_data,1);
    end
end

locations=[-9:1:9];

% fit via logistic regression
[a,b,c]=glmfit(suppDT.mag_dif,1-suppDT.choice,'binomial');

m = nanmean(av_ch,2);

% and fit the center via linear regression
[c,~,stats]  = glmfit(-4:4,m(6:14));

figure;hold on
plot(locations,nanmean(av_ch,2),'o','Color','k')

plot(locations,1./(1+exp(a(1)+a(2)*locations)),'--','Color','k')
plot([-4 4],[c(1) - 4* c(2),c(1) + 4*c(2)],'k')
xticks([-9:2:9]);
xlabel('magnitude difference (L - R), in ml') 
set(gca, 'FontSize',14)
ylabel('P(left)')
set(gcf,'color','w');
xticklabels(string([-4.5:4.5]))
box off


%% Fig S1C
% get the average timecourse for 1-10 drops in the orofacial cortex
[betas,names,regions] = fitTC(suppDTn,15,'d2');

% the intercept is the average (if there are no other predictors)
reg = '(Intercept)';
names{1}(1,1,1)='(Intercept)';

plottitle = '';
firstxlabel = "reward delivery";
regions = {'right orofacial cortex','left orofacial cortex'};

figure;hold on
for ii = 1:length(regions)
    subplot(2,2,ii);hold on
    plotTC(betas{ii},reg,names,plottitle,firstxlabel,suppDTn,15,"d2")

    title(regions{ii})
end

colormap(brewermap(10,'Blues'));
cbh = colorbar ;
cbh.Ticks = linspace(0.1, 1.9, 10) ;
cbh.TickLabels = num2cell(.5:.5:5);

colorTitleHandle = get(cbh,'Title');
titleString = 'juice (ml)';
set(colorTitleHandle ,'String',titleString);

% take the integral under the curve and plot it
monkey = [];
int1 = []; int2 = []; int3 = []; int4 = []; int5 = []; 
int6 = []; int7 = []; int8 = []; int9 = []; int10 = [];

t = -1:2.28/10:15;
search_window = [5:71];
for which_area = 1:2

    for m = 1:4
        for ii = 1:size(betas{which_area}.betasRew1,3)
            if ~isnan(betas{which_area}.betasRew1(1,m,ii,1))
                int1 = [int1;mean(betas{which_area}.betasRew1(search_window,m,ii,1))];
                int2 = [int2;mean(betas{which_area}.betasRew2(search_window,m,ii,1))];
                int3 = [int3;mean(betas{which_area}.betasRew3(search_window,m,ii,1))];
                int4 = [int4;mean(betas{which_area}.betasRew4(search_window,m,ii,1))];
                int5 = [int5;mean(betas{which_area}.betasRew5(search_window,m,ii,1))];
                int6 = [int6;mean(betas{which_area}.betasRew6(search_window,m,ii,1))];
                int7 = [int7;mean(betas{which_area}.betasRew7(search_window,m,ii,1))];
                int8 = [int8;mean(betas{which_area}.betasRew8(search_window,m,ii,1))];
                int9 = [int9;mean(betas{which_area}.betasRew9(search_window,m,ii,1))];
                int10 = [int10;mean(betas{which_area}.betasRew10(search_window,m,ii,1))];


                monkey = [monkey;m];
            end
        end
    end

    dat = table; dat.monkey = monkey;
    dat.int1 = int1; dat.int2 = int2; dat.int3 = int3; dat.int4 = int4;
    dat.int5 = int5; dat.int6 = int6; dat.int7 = int7; dat.int8 = int8;
    dat.int9 = int9; dat.int10 = int10;


    subplot(2,2,which_area+2); hold on;
    c = colormap(brewermap(10,'Blues'));

    all_means = nan(4,1);
    for ii = 1:10
        m_means = nan(4,1);
        for m = 1:4
            m_means(m) = mean(dat.(['int',num2str(ii)])(dat.monkey==m));
        end
        scatter(ii,mean(m_means),'MarkerEdgeColor','k','MarkerFaceColor',c(ii,:))
        all_means(ii) = mean(m_means);

    end
    
    % fit a line through the points
    [c,~,stats]  = glmfit(1:10,all_means);

    plot([1 10],[c(1) + c(2),c(1) + 10*c(2)],'k')

    xlabel('juice (ml)')
    ylabel('mean signal 0s to 15s (a.u.)')
    set(gca, 'FontSize',14)
    xticks(1:10)
    xticklabels(.5:.5:5)
end


%% Fig S1D
% get the average timecourse for 1-3 drops in the orofacial cortex
[betas,names,regions] = fitTC(DTn,15,'d1');

reg = '(Intercept)';
names{1}(1,1,1)='(Intercept)';

plottitle = '';
firstxlabel = "reward delivery";

regions = {'right orofacial cortex','left orofacial cortex'};

% plot it
figure;hold on
for ii = 1:length(regions)
    subplot(2,2,ii);hold on
    plotTC(betas{ii},reg,names,plottitle,firstxlabel,DT,15,"d2")

    title(regions{ii})
end

colormap(brewermap(10,'Blues'));
cbh = colorbar ;
cbh.Ticks = [0.5,1.1,1.7] ;
cbh.TickLabels = [1.5,3,4.5];

colorTitleHandle = get(cbh,'Title');
titleString = 'juice (ml)';
set(colorTitleHandle ,'String',titleString);

% also take the integral and plot it
monkey = []; int1 = []; int2 = []; int3 = []; 

t = -1:2.28/10:15;
search_window = [5:71];
for which_area = 1:2

    for m = 1:6
        for ii = 1:size(betas{which_area}.betasRew1,3)
            if ~isnan(betas{which_area}.betasRew1(1,m,ii,1))
                int1 = [int1;mean(betas{which_area}.betasRew1(search_window,m,ii,1))];
                int2 = [int2;mean(betas{which_area}.betasRew2(search_window,m,ii,1))];
                int3 = [int3;mean(betas{which_area}.betasRew3(search_window,m,ii,1))];


                monkey = [monkey;m];
            end
        end
    end

    dat = table;
    dat.monkey = monkey;
    dat.int1 = int1;
    dat.int2 = int2;
    dat.int3 = int3;


    subplot(2,2,which_area+2); hold on;
    c = colormap(brewermap(10,'Blues'));

    all_means = nan(3,1);
    for ii = 1:3
        m_means = nan(4,1);
        for m = 1:6
            m_means(m) = mean(dat.(['int',num2str(ii)])(dat.monkey==m));
        end
        scatter(ii,mean(m_means),'MarkerEdgeColor','k','MarkerFaceColor',c(ii*3,:))
        all_means(ii) = mean(m_means);

    end

    % fit line
    [c,~,stats]  = glmfit(1:3,all_means);

    plot([1 3],[c(1) + c(2),c(1) + 3*c(2)],'k')



    xlabel('juice (ml)')
    ylabel('mean signal 0s to 15s (a.u.)')
    set(gca, 'FontSize',14)
    xticks(1:3)
    xticklabels(1.5:1.5:4.5)
end



%% Fig S2
% run GLME1 for each session type and plot it

GLME1Changing = fitglme(DT(DT.Changing==1,:),'Lapse~VS+sRE+RRE+Position+TrialNo+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');
GLME1ChangingNoRRE = fitglme(DT(DT.Changing==1,:),'Lapse~sRE+VS+Position+TrialNo+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');
compare(GLME1ChangingNoRRE,GLME1Changing)

for Asub = 1:6
    GLME1ChangingMonkey{Asub} = fitglme(DT(DT.monkey==Asub&DT.Changing==1,:),'Lapse~VS+sRE+RRE+Position+TrialNo+(1|monkey:session)','distribution','binomial');
end

GLME1Stable = fitglme(DT(DT.Stable==1,:),'Lapse~VS+sRE+RRE+Position+TrialNo+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');
GLME1StableNoRRE = fitglme(DT(DT.Stable==1,:),'Lapse~sRE+VS+Position+TrialNo+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');
compare(GLME1StableNoRRE,GLME1Stable)

for Asub = 1:6
    GLME1StableMonkey{Asub} = fitglme(DT(DT.monkey==Asub&DT.Stable==1,:),'Lapse~VS+sRE+RRE+Position+TrialNo+(1|monkey:session)','distribution','binomial');
end

GLME1Equi = fitglme(DT(DT.Equiprobable==1,:),'Lapse~VS+sRE+RRE+Position+TrialNo+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');
GLME1EquiNoRRE = fitglme(DT(DT.Equiprobable==1,:),'Lapse~sRE+VS+Position+TrialNo+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','binomial','fitmethod','laplace');
compare(GLME1EquiNoRRE,GLME1Equi)

for Asub = 1:6
    GLME1EquiMonkey{Asub} = fitglme(DT(DT.monkey==Asub&DT.Equiprobable==1,:),'Lapse~VS+sRE+RRE+Position+TrialNo+(1|monkey:session)','distribution','binomial');
end

plotFitResultCoeff(GLME1Changing,GLME1Stable,GLME1Equi,GLME1ChangingMonkey,GLME1StableMonkey,GLME1EquiMonkey,4)

%% Fig S3A
% run response time GLMEs (no outlier trials)

GLME3 = fitglme(DT(DT.exclude==0,:),'RT~VS+sRE+RRE+Position+TrialNo+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME3nosRE = fitglme(DT(DT.exclude==0,:),'RT~VS+RRE+Position+TrialNo+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME3noVS = fitglme(DT(DT.exclude==0,:),'RT~sRE+RRE+Position+TrialNo+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME3noRRE = fitglme(DT(DT.exclude==0,:),'RT~sRE+VS+Position+TrialNo+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME3noTrialNo = fitglme(DT(DT.exclude==0,:),'RT~VS+sRE+RRE+Position+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');


GLME3monkey = {};
for Asub = 1:6
    GLME3monkey{Asub} = fitglme(DT(DT.monkey==Asub&DT.exclude==0,:),'RT~VS+sRE+RRE+Position+TrialNo+(1|monkey:session)','distribution','gamma','link','log');
end

plotFitResult(GLME3,GLME3monkey)

compare(GLME3nosRE,GLME3)
compare(GLME3noVS,GLME3)
compare(GLME3noRRE,GLME3)


%% Fig S3B
% run response time GLMEs (no outlier trials)

GLME4 = fitglme(DT(DT.exclude==0,:),'RT~VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME4not1 = fitglme(DT(DT.exclude==0,:),'RT~VS+t2+t3+t4+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME4not2 = fitglme(DT(DT.exclude==0,:),'RT~VS+t1+t3+t4+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME4not3 = fitglme(DT(DT.exclude==0,:),'RT~VS+t1+t2+t4+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME4not4 = fitglme(DT(DT.exclude==0,:),'RT~VS+t1+t2+t3+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME4not5 = fitglme(DT(DT.exclude==0,:),'RT~VS+t1+t2+t3+t4+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME4noVS = fitglme(DT(DT.exclude==0,:),'RT~t1+t2+t3+t4+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME4noTrialNo = fitglme(DT(DT.exclude==0,:),'RT~VS+t1+t2+t3+t4+t5+RRE+Position+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');


GLME4monkey = {};
for Asub = 1:6
    GLME4monkey{Asub} = fitglme(DT(DT.monkey==Asub&DT.exclude==0,:),'RT~VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo+(1|monkey:session)','distribution','gamma','link','log');
end

plotFitResult(GLME4,GLME4monkey)


%% Fig S3C
% run response time GLMEs (no outlier and repeat trials)

GLME5 = fitglme(DT(DT.exclude2==0,:),'RT~VS+sRE+RRE+Position+TrialNo+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME5nosRE = fitglme(DT(DT.exclude2==0,:),'RT~VS+RRE+Position+TrialNo+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME5noVS = fitglme(DT(DT.exclude2==0,:),'RT~sRE+RRE+Position+TrialNo+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME5noRRE = fitglme(DT(DT.exclude2==0,:),'RT~sRE+VS+Position+TrialNo+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME5noTrialNo = fitglme(DT(DT.exclude2==0,:),'RT~VS+sRE+RRE+Position+(VS+sRE+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');


GLME5monkey = {};
for Asub = 1:6
    GLME5monkey{Asub} = fitglme(DT(DT.monkey==Asub&DT.exclude2==0,:),'RT~VS+sRE+RRE+Position+TrialNo+(1|monkey:session)','distribution','gamma','link','log');
end
plotFitResult(GLME5,GLME5monkey)

compare(GLME5nosRE,GLME5)
compare(GLME5noVS,GLME5)
compare(GLME5noRRE,GLME5)

%% Fig S3D
% run response time GLMEs (no outlier and repeat trials)

GLME6 = fitglme(DT(DT.exclude2==0,:),'RT~VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME6not1 = fitglme(DT(DT.exclude2==0,:),'RT~VS+t2+t3+t4+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME6not2 = fitglme(DT(DT.exclude2==0,:),'RT~VS+t1+t3+t4+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME6not3 = fitglme(DT(DT.exclude2==0,:),'RT~VS+t1+t2+t4+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME6not4 = fitglme(DT(DT.exclude2==0,:),'RT~VS+t1+t2+t3+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME6not5 = fitglme(DT(DT.exclude2==0,:),'RT~VS+t1+t2+t3+t4+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME6noVS = fitglme(DT(DT.exclude2==0,:),'RT~t1+t2+t3+t4+t5+RRE+Position+TrialNo+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');
GLME6noTrialNo = fitglme(DT(DT.exclude2==0,:),'RT~VS+t1+t2+t3+t4+t5+RRE+Position+(VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo|monkey)+(1|monkey:session)','distribution','gamma','link','log','fitmethod','laplace');

save
GLME6monkey = {};
for Asub = 1:6
    GLME6monkey{Asub} = fitglme(DT(DT.monkey==Asub&DT.exclude2==0,:),'RT~VS+t1+t2+t3+t4+t5+RRE+Position+TrialNo+(1|monkey:session)','distribution','gamma','link','log');
end
plotFitResult(GLME6,GLME6monkey)


%% Fig S5A
% plot the z-stats within the ROI separately for each session type
plotROIzStatSessionType('PE_Cluster2',DTn,'right striatum',[-0.5 42.25 -inf inf])
[chi2 p] = getzStatp('PE_Cluster2',DTn)


%% Fig S5B
% plot the z-stats within the ROI separately for each session type
plotROIzStatSessionType('PE_Cluster1',DTn,'left striatum',[-0.5 42.25 -inf inf])
[chi2 p] = getzStatp('PE_Cluster1',DTn)


%% Fig S5C
% plot the z-stats within the ROI separately for each session type
plotROIzStatSessionType('SN_r',DTn,'right VTA/SN',[-0.5 42.25 -inf inf])
[chi2 p] = getzStatp('SN_r',DTn)


%% Fig S5D
% plot the z-stats within the ROI separately for each session type
plotROIzStatSessionType('SN_l',DTn,'left VTA/SN',[-0.5 42.25 -inf inf])
[chi2 p] = getzStatp('SN_l',DTn)


%% Fig S7
% run a regression on the timecourse in the VTA/SN with a range of learning
% rates
alpha_range = logspace(-3,0,30);
betas = cell(1,length(alpha_range));
names = cell(1,length(alpha_range));

parfor alpha = 1:length(alpha_range)
    % get the regressor with the new alpha
    DT_temp = arbitarayPE(alpha_range(alpha),DTn);
    % run the regression
    [betas{alpha},names{alpha},regions] = fitTC(DT_temp,10,'d3');
end

LR = 0.257; % the empirical alpha we fitted
t = -1:2.28/10:10;
search_window = [6:49];
tc = nan(6,length(t),length(alpha_range),2);
int = nan(71,length(alpha_range),2);
summary_m = nan(length(alpha_range),2);
summary_se = nan(length(alpha_range),2);

reg = 'sRPEtemp'; % what regressor to plot

% get the average area under the curve
for alpha = 1:length(alpha_range)
    r = 1;
    c = 1;
    monkey = [];
    for ID = 1:6
        for sess = 1:size(betas{alpha}{r}.betasRew,3)
            if ~isnan(mean(betas{alpha}{r}.betasRew(search_window,ID,sess,squeeze(names{1}{1}(1,1,:)==reg))))
            int(c,alpha,r) = mean(betas{alpha}{r}.betasRew(search_window,ID,sess,squeeze(names{1}{1}(1,1,:)==reg)));
            monkey = [monkey;ID];
            c = c + 1;
            end
        end
    end

    summary_m(alpha,r) = mean(int(:,alpha,r));
    summary_se(alpha,r) = std(int(:,alpha,r))/sqrt(71);
end

titles = {'VTA/SN'};
figure; hold on
% plot the average area under the curve for each learning rate
for r = 1
    shadedErrorBar(alpha_range,summary_m(:,r),summary_se(:,r))
    plot([LR LR],[0 max(summary_m(:,r))+0.005],'r')
    set(gca, 'XScale', 'log')
    title(titles(r))
    xlabel('alpha')
    ylabel('effect size (a.u.)')
    xticklabels(["0.001","0.01","0.1","1"])
end

ax = gca;
ax.FontSize = 16; 


%% Fig S8A
% plot the z-stats within the ROI separately for each session type
plotROIzStatSessionType('CurrentTwoDrop_Cluster1',DTn,'striatum',[-0.5 42.25 -inf inf])
[chi2 p] = getzStatp('CurrentTwoDrop_Cluster1',DTn)


%% Fig S8B
% plot the z-stats within the ROI separately for each session type
plotROIzStatSessionType('CurrentTwoDrop_lOFC',DTn,'plOFC',[-0.5 42.25 -inf inf])
[chi2 p] = getzStatp('CurrentTwoDrop_lOFC',DTn)


%% Fig S9
% plot the z-stats within the ROI separately for each session type
plotROIzStatSessionType('SSurprise_Cluster3',DTn,'lPFC',[-0.5 42.25 -inf inf])
[chi2 p] = getzStatp('SSurprise_Cluster3',DTn)


%% Fig S10C
% plot the z-stats within the ROI separately for each session type
plotROIzStatSessionType('ChangingEqui_striatum',DTn,'striatum',[-0.5 42.25 -inf inf])
% get the statistics for comparisons with the equiprobable session type
[chi2 p] = getzStatpContrast('ChangingEqui_striatum',DTn)


%% Fic S10D
% plot the z-stats within the ROI separately for each session type
plotROIzStatSessionType('ChangingEqui_lOFC',DTn,'plOFC',[-0.5 42.25 -inf inf])
% get the statistics for comparisons with the equiprobable session type
[chi2 p] = getzStatpContrast('ChangingEqui_lOFC',DTn)

