function edgeBoxesDemo2()
% Demo for Edge Boxes (please see readme.txt first).

%% load pre-trained edge detection model and set opts (see edgesDemo.m)
model=load('models/forest/modelBsds'); model=model.model;
model.opts.multiscale=1; model.opts.sharpen=2; model.opts.nThreads=4;

%% set up opts for edgeBoxes (see edgeBoxes.m)
opts = edgeBoxes;
opts.alpha = .65;     % step size of sliding window search
opts.beta  = .75;     % nms threshold for object proposals
opts.minScore = .01;  % min score of boxes to detect
opts.maxBoxes = 1e4;  % max number of boxes to detect

%% detect Edge Box bounding box proposals (see edgeBoxes.m)
filepath = '/data/images/everingham_IJCV_2010_pascalvoc/JPEGImages/2009_005302.jpg';
I = imread(filepath);
tic, [bbs,extra_info]=edgeBoxes2(I,model,opts); toc

visualizeBest(filepath, bbs, extra_info);
end

function visualizeBest(filepath, bbs, extra_info)
V = extra_info{1};
segIds = extra_info{2};
edgeArr = extra_info{3};
E = extra_info{4};
O = extra_info{5};
T = extra_info{6};

[maxf, maxidx] = getBestBB(filepath, bbs, segIds, edgeArr);

I = imread(filepath);
sz = size(I);
h = sz(1);  w = sz(2);

[ sp_seg ] = structedges_sp(I, E);
temp = rgb2gray(I);
% r = I(:,:,1);  g = I(:,:,2);  b = I(:,:,3);
r = temp; g = temp; b = temp;
avgclr = arrayfun(@(l) [mean(r(sp_seg == l)), mean(g(sp_seg == l)), mean(b(sp_seg == l))], 1:max(sp_seg(:)), 'UniformOutput', false);
avgclr = cell2mat(avgclr');
r = avgclr(sp_seg,1);  g = avgclr(sp_seg,2);  b = avgclr(sp_seg,3);
segclr = uint8(cat(3, reshape(r, [h,w]), reshape(g, [h,w]), reshape(b, [h,w])));

validIdsmask = segIds > 0;

c = uint8(winter(256) * 255);
v = linspace(0,1,256);
for idx = (maxidx(:))'
    fprintf('%d\n',idx);
    bbim = drawBox(bbs, idx, segIds, edgeArr);
    
    box = bbs(idx,:);

    r1 = clamp(box(2) + box(4), 1, h); 
    r0 = clamp(box(2), 1, h);
    c1 = clamp(box(1) + box(3), 1, w); 
    c0 = clamp(box(1), 1, w);
    
    o = bbim(validIdsmask);
    temp = bsxfun(@minus, o, v).^2;
    [~,clridx] = min(temp, [], 2);
    currclr = c(clridx,:);
    
    r = segclr(:,:,1);  g = segclr(:,:,2);  b = segclr(:,:,3);
    r(validIdsmask) = currclr(:,1);
    g(validIdsmask) = currclr(:,2);
    b(validIdsmask) = currclr(:,3);
    segclr = cat(3, r, g, b);
    
    bbim_sub = bbim(r0:r1, c0:c1);
    segclr_sub = segclr(r0:r1, c0:c1, :);
    
    imshow(segclr_sub, 'Parent',gca);
    pause;
end
end

function [ sp_seg ] = structedges_sp(orig_I, E)
    % set up opts for spDetect (see spDetect.m)
    opts = spDetect;
    opts.nThreads = 4;  % number of computation threads
    opts.k = 512;       % controls scale of superpixels (big k -> big sp)
    opts.alpha = .5;    % relative importance of regularity versus data terms
    opts.beta = .9;     % relative importance of edge versus color terms
    opts.merge = 0;     % set to small value to merge nearby superpixels at end

    [se_seg, V] = spDetect(orig_I, E, opts);

    sp_seg = fill_zeros(orig_I, se_seg);
    
    [sp_seg] = renumber_sps(sp_seg);
end

function [sp_seg] = renumber_sps(sp_seg) 
    % renumber superpixels, because it has missing indices
    sp_idxs = unique(sp_seg);
    mapping = zeros(sp_idxs(end),1);
    mapping(sp_idxs) = 1:length(sp_idxs);
    sp_seg = mapping(sp_seg);
end

function [sp_seg] = fill_zeros(orig_I, init_seg)
    addpath(fullfile(fileparts(which(mfilename)), '..', 'stein_boundaryprocessing'));
    
    % fill in the boundary pixels with the neighboring segment with
    % the closest RGB color
    sp_seg = fill_in_segmentation(orig_I, init_seg, 0, 4);
    
    % remove zeros (pixels which are still marked as boundaries)
    nghbrs = [-1 0 1];
    nghbrs = [nghbrs-size(sp_seg,1), nghbrs, nghbrs+size(sp_seg,1)];
    nghbrs(5) = [];
    zero_locs = find(sp_seg == 0);
    nghbrs = bsxfun(@plus, zero_locs, nghbrs);
    nghbrs = sp_seg(nghbrs);
    % setting it to mode, which is the most frequently occuring
    % superpixel in the neighborhood
    sp_seg(zero_locs) = uint16(mode(single(nghbrs), 2));
end

function [v] = clamp(v, ll, ul)
if v < ll
    v = ll;
elseif v > ul;
    v = ul;
end
end