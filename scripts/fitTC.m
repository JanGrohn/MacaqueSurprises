function [betas,names,regions] = fitTC(DT,window,dat)

% what regions to look at
if dat == "d1"
    regions = {'peak5_al','peak6_al','PE_Cluster2','PE_Cluster1','SN_l','SN_r'};
elseif dat == "d2"
   regions = {'peak5_al','peak6_al'}; 
elseif dat == "d3"
    regions = {{'SN_l','SN_r'}};
end


betas = {};
names = {};

n = length(regions);

for region = 1:n
    display(['fitting region ',num2str(region),'/',num2str(n)]);
    % what analysis to run
    if dat == "d1"
        [betas{region},names{region}] = fitTCIndSess(DT,regions{region},window);
    elseif dat == "d2"
        [betas{region},names{region}] = fitTCIndSessSupp(DT,regions{region},window);
    elseif dat == "d3"
        [betas{region},names{region}] = fitTCIndSessAlpha(DT,regions{region},window);
    end

end

display("end")

end
