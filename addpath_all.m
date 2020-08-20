% adds all folders and subfolders to MATLAB path
files = dir;
folders = files([files.isdir]);
for ii = 3:size(folders,1)
    addpath(genpath(folders(ii).name));
end
clear files folders ii