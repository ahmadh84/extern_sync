% This demo shows how to use the software described in our IJCV paper: 
%   Selective Search for Object Recognition,
%   J.R.R. Uijlings, K.E.A. van de Sande, T. Gevers, A.W.M. Smeulders, IJCV 2013
%%
addpath('Dependencies');

fprintf('Demo of how to run the code for:\n');
fprintf('   J. Uijlings, K. van de Sande, T. Gevers, A. Smeulders\n');
fprintf('   Segmentation as Selective Search for Object Recognition\n');
fprintf('   IJCV 2013\n\n');

% Compile anisotropic gaussian filter
if(~exist('anigauss'))
    fprintf('Compiling the anisotropic gauss filtering of:\n');
    fprintf('   J. Geusebroek, A. Smeulders, and J. van de Weijer\n');
    fprintf('   Fast anisotropic gauss filtering\n');
    fprintf('   IEEE Transactions on Image Processing, 2003\n');
    fprintf('Source code/Project page:\n');
    fprintf('   http://staff.science.uva.nl/~mark/downloads.html#anigauss\n\n');
    mex Dependencies/anigaussm/anigauss_mex.c Dependencies/anigaussm/anigauss.c -output anigauss
end

if(~exist('mexCountWordsIndex'))
    mex Dependencies/mexCountWordsIndex.cpp
end

% Compile the code of Felzenszwalb and Huttenlocher, IJCV 2004.
if(~exist('mexFelzenSegmentIndex'))
    fprintf('Compiling the segmentation algorithm of:\n');
    fprintf('   P. Felzenszwalb and D. Huttenlocher\n');
    fprintf('   Efficient Graph-Based Image Segmentation\n');
    fprintf('   International Journal of Computer Vision, 2004\n');
    fprintf('Source code/Project page:\n');
    fprintf('   http://www.cs.brown.edu/~pff/segment/\n');
    fprintf('Note: A small Matlab wrapper was made.\n');
%     fprintf('   
    mex Dependencies/FelzenSegment/mexFelzenSegmentIndex.cpp -output mexFelzenSegmentIndex;
end

%%
% Parameters. Note that this controls the number of hierarchical
% segmentations which are combined.
colorTypes = {'Hsv', 'Lab', 'RGI', 'H', 'Intensity'};
% colorType = colorTypes{1}; % Single color space for demo

% Here you specify which similarity functions to use in merging
simFunctionHandles = {@SSSimColourTextureSizeFillOrig, @SSSimTextureSizeFill, @SSSimBoxFillOrig, @SSSimSize};
% simFunctionHandles = simFunctionHandles(1:2); % Two different merging strategies

% Thresholds for the Felzenszwalb and Huttenlocher segmentation algorithm.
% Note that by default, we set minSize = k, and sigma = 0.8.
% k = 200; % controls size of segments of initial segmentation. 
% minSize = k;
ks = [50 100 150 300]; % controls size of segments of initial segmentation. 
sigma = 0.8;

% After segmentation, filter out boxes which have a width/height smaller
% than minBoxWidth (default = 20 pixels).
minBoxWidth = 20;

% As an example, use a single image
images = {'/home/ahumayun/Dropbox/NFoundSeg/JPEGImages/2008_007998.jpg'};
im = imread(images{1});

% % Perform Selective Search
% [boxes blobIndIm blobBoxes hierarchy] = Image2HierarchicalGrouping(im, sigma, k, minSize, colorType, simFunctionHandles);
% boxes = BoxRemoveDuplicates(boxes);
% 
% % Show boxes
% % ShowRectsWithinImage(boxes, 5, 5, im);
% 
% % Show blobs which result from first similarity function
% hBlobs = RecreateBlobHierarchyIndIm(blobIndIm, blobBoxes, hierarchy{1});
% ShowBlobs(hBlobs, 5, 5, im);

totalTime = 0;
idx = 1;
for j=1:length(ks)
    k = ks(j); % Segmentation threshold k
    minSize = k; % We set minSize = k
    for n = 1:length(colorTypes)
        colorType = colorTypes{n};
        tic;
        [boxesT{idx} blobIndIm blobBoxes hierarchy priorityT{idx}] = Image2HierarchicalGrouping(im, sigma, k, minSize, colorType, simFunctionHandles);
        totalTime = totalTime + toc;
        idx = idx + 1;
    end
end
boxes = cat(1, boxesT{:}); % Concatenate boxes from all hierarchies
priority = cat(1, priorityT{:}); % Concatenate priorities

% Do pseudo random sorting as in paper
priority = priority .* rand(size(priority));
[priority, sortIds] = sort(priority, 'ascend');
boxes = boxes(sortIds,:);
fprintf('\n');

tic
boxes = FilterBoxesWidth(boxes, minBoxWidth);
boxes = BoxRemoveDuplicates(boxes);
totalTime = totalTime + toc;

fprintf('Time for this image %.2f ...\n', totalTime);
