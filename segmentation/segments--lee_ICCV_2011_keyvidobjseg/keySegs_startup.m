function keySegs_startup
%KEYSEGS_STARTUP Summary of this function goes here
%   Detailed explanation goes here

    script_root = fileparts(which(mfilename));
    addpath(genpath(fullfile(script_root, '/code/segmentation/')));
    addpath(genpath(fullfile(script_root, '/code/getHypotheses/')));
    addpath(genpath(fullfile(script_root, '/code/util/')));
    addpath(fullfile(script_root, '..', '..', 'metrics', 'hua_pwmetric'));
    addpath(fullfile(script_root, '..', ...
        'merges--malisiewicz_BMVC_2007_slidingsegments'));
    addpath(fullfile(script_root, '..', '..', 'machine_learning', ...
        'calinon_TrscSMC_2007_GMMGMR-v2.0'));
    addpath(fullfile(script_root, '..', '..', 'optimization', ...
        'maxflow--boykov_PAMI_2001_GCmex1.5'));
end

