% hypotheses_demo
%
% do run the demo, download this data first: 
% vision.cs.utexas.edu/projects/keysegments/data.tgz
%
clear;

script_root = fileparts(which(mfilename));
addpath(genpath(fullfile(script_root, '/code/getHypotheses/')));
addpath(genpath(fullfile(script_root, '/code/util/')));
addpath(fullfile(script_root, '..', '..', 'metrics', 'hua_pwmetric'));
addpath(fullfile(script_root, '..', 'merges--malisiewicz_BMVC_2007_slidingsegments'));

%%%%%%%%%%
skip = 1;
numTopRegions = 10;
knn = 5;

IMSHOW = 1;
%%%%%%%%%%

videoName = 'girl'; hypothesisNum = 1;
% videoName = 'parachute'; hypothesisNum = 1;
% videoName = 'birdfall2'; hypothesisNum = 1;
% videoName = 'monkeydog'; hypothesisNum = 3;
% videoName = 'cheetah'; hypothesisNum = 2;
% videoName = 'penguin'; hypothesisNum = 1;

imdir = ['~/data/Images/SegTrack_201102/' videoName '/'];
regiondir = ['data/regionProposals/' videoName '/'];
opticalflowdir = ['data/opticalFlow/' videoName '/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get hypotheses (i.e., keysegment clusters)
tic;
[A Segs imNdx scores hypotheses] = getRegionAffinitiesMotion(imdir, regiondir, opticalflowdir, skip, numTopRegions, knn);
selind = hypotheses(hypothesisNum).selind;
toc;

% display the key-segments in hypothesis = hypothesisNum
if IMSHOW
    figure(1); clf;
    for i=1:length(selind)
        im = imread([imdir Segs(selind(i)).imname]);
        I2 = segImage2(im2double(im),double(Segs(selind(i)).proposal),2);
        imagesc(I2);
        pause;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

