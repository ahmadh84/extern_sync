function computeRegionProposals(videoName, data_dir, num_tries)

    curr_dir = fileparts(which(mfilename));

    videoDir = [fullfile(data_dir, videoName) '/'];  % path to the video
    
    addpath(fullfile(curr_dir, '..', '..', '..', 'segments--endres_ECCV_2010_objsegproposals'));

    outDir = [curr_dir, '/../../data/regionProposals/' videoName '/'];    % output directory
    if exist('num_tries','var') ~= 1
        num_tries = 1;
    end
    
   for idx = 1:num_tries
       curr_outDir = outDir;
       if num_tries > 1
           curr_outDir = [curr_outDir, sprintf('test_%02d', idx), '/'];
       end
       runRegionProposals(videoDir, curr_outDir);
   end
end

function runRegionProposals(videoDir, outDir)
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
        save('-v7',[outDir  imname '.mat'], 'proposals', 'superpixels', 'image_data', 'unary');
    end
end