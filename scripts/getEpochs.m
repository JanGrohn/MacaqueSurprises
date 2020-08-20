function data = getEpochs(experiment,region, data, window)

%% Upsample timecourse 
% Presets for upsampling
TR=2.28;  % from scan
upsample=10; % upsampling factor


%% Load and upsample the timecourse file

rewEpochs = [];

for Asub = 1:max(data.monkey)
    sessions = unique(data.session(data.monkey==Asub))';
    for sess = sessions

        name = data.name(data.monkey==Asub&data.session==sess);
        sessID = data.sessID(data.monkey==Asub&data.session==sess);
        loadpathtimc = fullfile('timeseries',name(1),sessID(1),strcat(region,'.txt'));
        if strcmp(experiment,'d1')
            nuisanceName = 'nuisance';
        elseif strcmp(experiment,'d2') 
            nuisanceName = 'confound_EVs_2_5std_noMEL';
        end
        loadpathconf = fullfile('timeseries',name(1),sessID(1),[nuisanceName,'.txt']);
        timecourse=load(loadpathtimc); % load the timecourse (BOLD responses at the measured times, each measurement is TR apart from one before)
        confounds=dlmread(loadpathconf); % get nuisance regressors
        
        timecourse=zscore(timecourse); % normalise
        
        [~,~,stats]=glmfit(confounds,timecourse); % regress out nuisance
        timecourse = stats.resid;
        
        timecourse=zscore(timecourse); % normalise again

        t=0:length(timecourse)-1;
        t_ups=0:(1/upsample):length(timecourse)-1;
        timecourse=spline(t,timecourse,t_ups); % upsample the timecourse
        
        
        % reward period
        rewEpochs = [rewEpochs;cutEpochs(data.rewOnset(data.monkey==Asub&data.session==sess),timecourse)];

    end
end

data.rewEpochs = rewEpochs;

function out = cutEpochs(onsets,timecourse)
    
trind=int32(ceil((onsets/TR)*upsample)-(upsample/2)); 
tps=-1:TR/upsample:window;

out =[];


for trial=1:length(trind)
    if trind(trial)-length(-1:TR/upsample:0)<=0
        out=[out; [nan(1,abs(trind(trial)-length(-1:TR/upsample:0))),timecourse(1:length(tps)-abs(trind(trial)-length(-1:TR/upsample:0)))]];
    elseif (trind(trial)+length(0:TR/upsample:window)-1) > length(timecourse)
        out=[out; [timecourse(trind(trial)-length(-1:TR/upsample:0):length(timecourse)),nan(1,length(tps)-length(trind(trial)-length(-1:TR/upsample:0):length(timecourse)))]];
    else
        out=[out; timecourse(trind(trial)-length(-1:TR/upsample:0):(trind(trial)+length(0:TR/upsample:window)-1))];
    end
end

% normalise
out=out-nanmean(out(:));
out=out./nanstd(out(:));

end
end