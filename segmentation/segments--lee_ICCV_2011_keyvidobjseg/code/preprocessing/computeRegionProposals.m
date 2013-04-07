function computeRegionProposals(videoName, data_dir)

    curr_dir = fileparts(which(mfilename));

    addpath(fullfile(curr_dir, '..', '..', '..', 'segments--endres_ECCV_2010_objsegproposals'));

    outDir = [curr_dir, '/../../data/regionProposals/' videoName '/'];    % output directory
    videoDir = [fullfile(data_dir, videoName) '/'];  % path to the video

    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end

    d = dir([videoDir '/*.jpg']);
    if isempty(d)
        d = dir([videoDir '/*.png']);
    end
    if isempty(d)
        d = dir([videoDir '/*.bmp']);
    end

    for i=1:length(d)
        imname = d(i).name;

        [proposals superpixels image_data unary] = generate_proposals([videoDir imname]);
        save('-v7',[outDir  imname '.mat'], 'proposals', 'superpixels', 'unary');
    end
end