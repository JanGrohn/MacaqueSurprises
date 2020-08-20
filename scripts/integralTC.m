function [stats] = integralTC(betas,which_area)

monkey = [];
int1 = [];
int2 = [];
int3 = [];
t = -1:2.28/10:10;
search_window = [6:49];
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

fitInt1 = fitlme(dat,'int1~1+(1|monkey)');
fitInt2 = fitlme(dat,'int2~1+(1|monkey)');
fitInt3 = fitlme(dat,'int3~1+(1|monkey)');

fitInt1n = fitlme(dat,'int1~-1+(1|monkey)');
fitInt2n = fitlme(dat,'int2~-1+(1|monkey)');
fitInt3n = fitlme(dat,'int3~-1+(1|monkey)');

stats1 = compare(fitInt1n,fitInt1);
stats2 = compare(fitInt2n,fitInt2);
stats3 = compare(fitInt3n,fitInt3);

stats = nan(3,2);
stats(1,1) = stats1.LRStat(2);
stats(2,1) = stats2.LRStat(2);
stats(3,1) = stats3.LRStat(2);
stats(1,2) = stats1.pValue(2);
stats(2,2) = stats2.pValue(2);
stats(3,2) = stats3.pValue(2);
end