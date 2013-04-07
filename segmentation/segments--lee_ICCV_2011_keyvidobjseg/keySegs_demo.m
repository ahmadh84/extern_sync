function keySegs_demo(videoName, data_dir, precompute_dir, hypothesisNum)
% keySegs_demo
%
% do run the demo, download this data first: 
% vision.cs.utexas.edu/projects/keysegments/data.tgz
%

    script_root = fileparts(which(mfilename));
    addpath(genpath(fullfile(script_root, '/code/segmentation/')));
    addpath(genpath(fullfile(script_root, '/code/getHypotheses/')));
    addpath(genpath(fullfile(script_root, '/code/util/')));
    addpath(fullfile(script_root, '..', '..', 'metrics', 'hua_pwmetric'));
    addpath(fullfile(script_root, '..', 'merges--malisiewicz_BMVC_2007_slidingsegments'));
    addpath(fullfile(script_root, '..', '..', 'machine_learning', 'calinon_TrscSMC_2007_GMMGMR-v2.0'));
    addpath(fullfile(script_root, '..', '..', 'optimization', 'maxflow--boykov_PAMI_2001_GCmex1.5'));
    
    %% parameters
    %%%%%%%%%%
    BPLR_ = 1;

    param.fgNbStates = 5;
    param.bgNbStates = 5;
    param.alpha_ = 0.5;
    if BPLR_ == 1
        param.gamma_ = 4;
    else
        param.gamma_ = 50;
    end

    param.skip = 1;
    param.numTopRegions = 10;
    param.knn = 5;

    param.hypothesisNum = hypothesisNum;
    
    param.imdir = [fullfile(data_dir, videoName) '/'];
    param.opticalflowdir = [fullfile(precompute_dir, 'opticalFlow', videoName) '/'];
    param.regiondir = [fullfile(precompute_dir, 'regionProposals', videoName) '/'];
    param.bplrdir = [fullfile(precompute_dir, 'bplr', videoName) '/'];
    %%%%%%%%%%

    [our_mask, Segs, selind] = videoSegmentation(BPLR_, param);
    % save('-v7',['./test/' videoName '.mat'],'our_mask', 'param', 'Segs', 'selind');
    % save('-v7',['./test/' videoName '_colorOnly.mat'],'our_mask', 'param', 'Segs', 'selind');

    evaluateSegTrackResults(param.imdir,our_mask,Segs,selind);
end