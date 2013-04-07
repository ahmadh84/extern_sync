function computeOpticalFlow(videoName, data_dir)
    curr_dir = fileparts(which(mfilename));

    addpath(fullfile(curr_dir, '..', '..', '..', '..', 'temporal_correspondences', 'flow--brox_PAMI_2010_ldof'));

    outDir = [curr_dir, '/../../data/opticalFlow/' videoName '/'];  % output directory
    % videoDir = ['~/data/Images/SegTrack_201102/' videoName '/'];    % path to the video
    videoDir = [fullfile(data_dir, videoName) '/'];    % path to the video

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

    for i=1:length(d)-1

        imname1 = d(i).name;
        imname2 = d(i+1).name;

        im1 = double(imread([videoDir imname1]));
        im2 = double(imread([videoDir imname2]));

        tic;
        flow = mex_LDOF(im1,im2);
        toc;

    %     figure(1); clf;
    %     subplot(2,2,1); imshow([videoDir imname1]);
    %     subplot(2,2,2); imshow([videoDir imname2]);
    %     subplot(2,2,3); imagesc(flow(:,:,1));
    %     subplot(2,2,4); imagesc(flow(:,:,2));
    %     pause;

        vx = flow(:,:,1);
        vy = flow(:,:,2);

        outFile = [outDir imname1 '_to_' imname2 '.opticalflow.mat'];
        save('-v7',outFile,'vx','vy');

        fprintf(1, 'LDOF computed (%d/%d) = %s\n', i, length(d)-1, outFile);
    end
end
