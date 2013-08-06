function keySegs_demo(videoName, data_dir, precompute_dir, hypothesisNum)
% keySegs_demo
%
% do run the demo, download this data first: 
% vision.cs.utexas.edu/projects/keysegments/data.tgz
%

    keySegs_startup;
    
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