% Script to generate selective search results in a way that its easy to get
% the segmentation proposals from the results. Plus it can be instructed to
% generate regions according to the three settings prescribed in their IJCV
% paper: single, fast, and quality.
%
% @authors:     Ahmad Humayun
% @contact:     ahumayun@cc.gatech.edu
% @affiliation: Georgia Institute of Technology
% @date:        Fall 2014

% This demo shows how to use the software described in our IJCV paper: 
%   Selective Search for Object Recognition,
%   J.R.R. Uijlings, K.E.A. van de Sande, T. Gevers, A.W.M. Smeulders, IJCV 2013
%%
function [segs, seg_obj] = selective_search_segs(im, params)
    root_dir = fileparts(which(mfilename));
    addpath(fullfile(root_dir, 'Dependencies'));
    
    if exist('params','var')==1 && ischar(params) 
        fprintf('Running selective search with ''%s'' setting\n', ...
                params);
    end
    
    if exist('params','var')==0 || (ischar(params) && strcmp(params,'single'))
        params = single_params();
    elseif ischar(params) && strcmp(params,'fast')
        params = fast_params();
    elseif ischar(params) && strcmp(params,'quality')
        params = quality_params();
    else
        error('Invalid params to Selective Search');
    end
    
    % As an example, use a single image
    % images = {'/home/ahumayun/Dropbox/NFoundSeg/JPEGImages/2008_007998.jpg'};
    % im = imread(images{1});
    
    seg_obj.segm_params = params;
    
    num_cfgs = length(params.colorTypes)*numel(params.ks);
    blobIndIm = zeros(size(im,1), size(im,2), num_cfgs);
    fprintf('Number of strategies: %d\n', num_cfgs*length(params.simFunctionHandles));
    
    start_time = tic();
    hiergrptime = 0;
    idx = 1;
    for j=1:length(params.ks)
        k = params.ks(j); % Segmentation threshold k
        minSize = k; % We set minSize = k
        for n = 1:length(params.colorTypes)
            colorType = params.colorTypes{n};
            t = tic;
            [boxesT{idx}, blobIndIm(:,:,idx), blobBoxes{idx}, hierarchies{idx}, ...
             priorityT{idx}, simfunc_idxs{idx}, order_idxs{idx}] = ...
                Image2HierarchicalGrouping(im, params.sigma, k, ...
                                           minSize, colorType, ...
                                           params.simFunctionHandles);
            hiergrptime = hiergrptime + toc(t);
            idx = idx + 1;
        end
    end
    seg_obj.timings.hier_grp_time = hiergrptime;
    
    num_segs = cellfun(@numel, priorityT);
    cfg_idxs = duplicateElems(1:num_cfgs, num_segs);

    bb_boxes = cat(1, boxesT{:}); % Concatenate boxes from all hierarchies
    priority = cat(1, priorityT{:}); % Concatenate priorities
    simfunc_idxs = cat(1, simfunc_idxs{:}); % concatenate function similarity idxs
    
    offset_idx = [0, cumsum(num_segs(1:end-1))];
    order_idxs = arrayfun(@(i) order_idxs{i}+offset_idx(i), 1:numel(order_idxs), 'UniformOutput',false);
    order_idxs = cat(1, order_idxs{:});
    
    seg_obj.num_segs.init_segs = num_segs;
    
    % Do pseudo random sorting as in paper
    t = tic;
    priority = priority .* rand(size(priority));
    [priority, sortIds] = sort(priority, 'ascend');
    bb_boxes = bb_boxes(sortIds,:);
    cfg_idxs = cfg_idxs(sortIds);
    simfunc_idxs = simfunc_idxs(sortIds);
    order_idxs = order_idxs(sortIds);
    seg_obj.timings.random_sort_time = toc(t);
    
    t = tic;
    [bb_boxes, idsGood] = FilterBoxesWidth(bb_boxes, params.minBoxWidth);
    priority = priority(idsGood);
    cfg_idxs = cfg_idxs(idsGood);
    simfunc_idxs = simfunc_idxs(idsGood);
    order_idxs = order_idxs(idsGood);
    seg_obj.num_segs.after_boxwidth_filter = accumarray(cfg_idxs(:), 1)';
    seg_obj.timings.boxwidth_filter_time = toc(t);
    
    t = tic;
    [bb_boxes, uniqueIdx] = BoxRemoveDuplicates(bb_boxes);
    priority = priority(uniqueIdx);
    cfg_idxs = cfg_idxs(uniqueIdx);
    simfunc_idxs = simfunc_idxs(uniqueIdx);
    order_idxs = order_idxs(uniqueIdx);
    seg_obj.num_segs.after_repeat_remove = accumarray(cfg_idxs(:), 1)';
    seg_obj.timings.box_similar_filter_time = toc(t);

    seg_obj.num_segs.final_num_segs = seg_obj.num_segs.after_repeat_remove;
    
    % store data which can be used to generate segments
    segs.scores = priority;
    segs.cfg_idxs = cfg_idxs;
    segs.simfunc_idxs = simfunc_idxs;
    segs.order_idxs = order_idxs;
    segs.sp_maps = blobIndIm;
    segs.hierarchies = hierarchies;
    segs.bb_boxes = bb_boxes;
    segs.sp_blob_boxes = blobBoxes;
    
    seg_obj.timings.total_seg_time = toc(start_time);
    
    fprintf('Time for this image %.2f ...\n', ...
            seg_obj.timings.total_seg_time);
end


function [params] = single_params()
    params = quality_params();
    
    params.colorTypes = params.colorTypes(1);   % HSV
    params.simFunctionHandles = params.simFunctionHandles(1);  % C+T+S+F
    params.ks = 100; % controls size of segments of initial segmentation. 
end

function [params] = fast_params()
    params = quality_params();
    
    params.colorTypes = params.colorTypes(1:2); % HSV and Lab
    params.simFunctionHandles = params.simFunctionHandles(1:2); % C+T+S+F,  T+S+F
    params.ks = params.ks(1:2);
end

function [params] = quality_params()
    %%
    % Parameters. Note that this controls the number of hierarchical
    % segmentations which are combined.
    params.colorTypes = {'Hsv', 'Lab', 'RGI', 'H', 'Intensity'};

    % Here you specify which similarity functions to use in merging
    %                                                            T+S+F
    params.simFunctionHandles = {@SSSimColourTextureSizeFillOrig, ... % C+T+S+F
                                 @SSSimTextureSizeFill, ...           % T+S+F
                                 @SSSimBoxFillOrig, ...               % F
                                 @SSSimSize};                         % S

    % Thresholds for the Felzenszwalb and Huttenlocher segmentation algorithm.
    % Note that by default, we set minSize = k, and sigma = 0.8.
    % k = 200; % controls size of segments of initial segmentation. 
    % minSize = k;
    params.ks = [50 100 150 300]; % controls size of segments of initial segmentation. 
    params.sigma = 0.8;

    % After segmentation, filter out boxes which have a width/height smaller
    % than minBoxWidth (default = 20 pixels).
    params.minBoxWidth = 20;
end
