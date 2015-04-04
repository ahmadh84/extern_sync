% ------------------------------------------------------------------------ 
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
%  University of California Berkeley (UCB) - USA
% 
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  Pablo Arbelaez <arbelaez@berkeley.edu>
%  June 2014
% ------------------------------------------------------------------------ 
% This file is part of the MCG package presented in:
%    Arbelaez P, Pont-Tuset J, Barron J, Marques F, Malik J,
%    "Multiscale Combinatorial Grouping,"
%    Computer Vision and Pattern Recognition (CVPR) 2014.
% Please consider citing the paper if you use this code.
% ------------------------------------------------------------------------
% This function computes the MCG candidates given an image.
%  INPUT:
%  - image : Input image
%  - mode  : It can be: + 'fast'     (SCG in the paper)
%                       + 'accurate' (MCG in the paper)
%  - compute_masks : Compute the candidates maks [0] or not [1]. Note that
%                    it is very time consuming. Otherwise use the labels 
%                    as shown in the demo.
%
%  OUTPUT:
%  - candidates : Struct containing the following fields
%          + superpixels : Label matrix of the superpixel partition
%          + labels : Cell containing the superpixel labels that form
%                     each of the candidates
%          + scores : Score of each of the ranked candidates
%          + masks  : 3D boolean matrix containing the masks of the candidates
%                     (only if compute_masks==1)
%  - ucm2       : Ultrametric Contour Map from which the candidates are
%                 extracted
%  - bboxes     : Bounding boxes of the candidates (up,left,down,right)
%                 See 'bboxes' folder for functions to work with them
%
%  DEMO:
%  - See demos/demo_im2mcg.m
% ------------------------------------------------------------------------
function [candidates, ucm2, seg_obj] = im2mcg(image,mode,compute_masks,preload)
if nargin<2
    mode = 'fast';
end
if nargin<3
    compute_masks = 0;
end

if exist('preload','var') ~= 1
    preload = struct;
end

% Load pre-trained Structured Forest model
if ~isfield(preload, 'sf_model')
    preload.sf_model = ...
        loadvar(fullfile(mcg_root, 'datasets', 'models', 'sf_modelFinal.mat'),'model');
end

start_time = tic();
total_load_time = 0;

% Level of overlap to erase duplicates
J_th = 0.95;

% Max margin parameter
theta = 0.7;

t = tic();
if strcmp(mode,'fast')
    % Which scales to work on (MCG is [2, 1, 0.5], SCG is just [1])
    scales = 1;

    % Get the hierarchies at each scale and the global hierarchy
    [ucm2] = img2ucms(image, preload.sf_model, scales);
    all_ucms = ucm2;
    
    t2 = tic();
    pareto_n_cands_fp = fullfile(mcg_root, 'datasets', 'models', 'scg_pareto_point_train2012.mat');
    % Load pre-trained pareto point
    if ~isfield(preload, 'pareto_n_cands')
        preload.pareto_n_cands = loadvar(pareto_n_cands_fp, 'n_cands');
    end

    rf_regressor_fp = fullfile(mcg_root, 'datasets', 'models', 'scg_rand_forest_train2012.mat');
    % Load pre-trained random forest regresssor for the ranking of candidates
    if ~isfield(preload, 'rf_regressor')
        preload.rf_regressor = loadvar(rf_regressor_fp,'rf');
    end
    total_load_time = toc(t2);
    
elseif strcmp(mode,'accurate')
    % Which scales to work on (MCG is [2, 1, 0.5], SCG is just [1])
    scales = [2, 1, 0.5];

    % Get the hierarchies at each scale and the global hierarchy
    [ucm2,ucms,times] = img2ucms(image, preload.sf_model, scales);
    all_ucms = cat(3,ucm2,ucms(:,:,3),ucms(:,:,2),ucms(:,:,1)); % Multi, 0.5, 1, 2

    t2 = tic();
    pareto_n_cands_fp = fullfile(mcg_root, 'datasets', 'models', 'mcg_pareto_point_train2012.mat');
    % Load pre-trained pareto point
    if ~isfield(preload, 'pareto_n_cands')
        preload.pareto_n_cands = loadvar(pareto_n_cands_fp, 'n_cands');
    end

    rf_regressor_fp = fullfile(mcg_root, 'datasets', 'models', 'mcg_rand_forest_train2012.mat');
    % Load pre-trained random forest regresssor for the ranking of candidates
    if ~isfield(preload, 'rf_regressor')
        preload.rf_regressor = loadvar(rf_regressor_fp, 'rf');
    end
    total_load_time = toc(t2);
else
    error('Unknown mode for MCG: Possibilities are ''fast'' or ''accurate''')
end
seg_obj.timings.ucm_time = toc(t) - total_load_time;
% ------------------------------------

% Transform ucms to hierarchies (dendogram) and put them all together
t = tic();
n_hiers = size(all_ucms,3);
lps = [];
ms  = cell(n_hiers,1);
ths = cell(n_hiers,1);
for ii=1:n_hiers
    % Transform the UCM to a hierarchy
    curr_hier = ucm2hier(all_ucms(:,:,ii));
    ths{ii}.start_ths = curr_hier.start_ths';
    ths{ii}.end_ths   = curr_hier.end_ths';
    ms{ii}            = curr_hier.ms_matrix;
    lps = cat(3, lps, curr_hier.leaves_part);
end
% detect empty hierarchies
empty_hier = cellfun(@isempty, ms);
seg_obj.timings.ucm_to_hier_time = toc(t);

% Get full cands, represented on a fused hierarchy
t = tic();
[f_lp,f_ms,cands,start_ths,end_ths] = full_cands_from_hiers(lps,ms,ths,preload.pareto_n_cands);
seg_obj.timings.cand_from_hier_time = toc(t);
seg_obj.num_segs.init_segs = size(cands,1);

% Hole filling and complementary candidates
t = tic();
if ~isempty(f_ms)
    [cands_hf, cands_comp] = hole_filling(double(f_lp), double(f_ms), cands); %#ok<NASGU>
else
    cands_hf = cands;
    cands_comp = cands; %#ok<NASGU>
end
seg_obj.timings.hole_filling_time = toc(t);
seg_obj.num_segs.hole_filled_segs = size(cands_hf,1);

% Select which candidates to keep (Uncomment just one line)
cands = cands_hf;                       % Just the candidates with holes filled
% cands = [cands_hf; cands_comp];         % Holes filled and the complementary
% cands = [cands; cands_hf; cands_comp];  % All of them
        
% Compute base features
t = tic();
b_feats = compute_base_features(f_lp, f_ms, all_ucms);
b_feats.start_ths = start_ths;
b_feats.end_ths   = end_ths;
b_feats.im_size   = size(f_lp);
seg_obj.timings.filter_feats_comp_time = toc(t);

% Filter by overlap
t = tic();
red_cands = mex_fast_reduction(cands-1,b_feats.areas,b_feats.intersections,J_th);
seg_obj.timings.seg_similar_filter_time = toc(t);
seg_obj.num_segs.after_repeat_remove = size(red_cands,1);

% Compute full features on reduced cands
t = tic();
[feats, bboxes] = compute_full_features(red_cands,b_feats);
seg_obj.timings.rank_feats_time = toc(t);

% Rank candidates
t = tic();
class_scores = regRF_predict(feats,preload.rf_regressor);
[scores, ids] = sort(class_scores,'descend');
red_cands = red_cands(ids,:);
bboxes = bboxes(ids,:);
if isrow(scores)
    scores = scores';
end
seg_obj.timings.rank_time = toc(t);

% Max margin
t = tic();
[new_ids, mm_scores] = mex_max_margin(red_cands-1,scores,b_feats.intersections,theta); %#ok<NASGU>
cand_labels = red_cands(new_ids,:);
candidates.scores = scores(new_ids);
bboxes = bboxes(new_ids,:); 
seg_obj.timings.max_margin_time = toc(t);

% Filter boxes by overlap
J_th_mex_box_red = 0.95;
red_bboxes = mex_box_reduction(bboxes, J_th_mex_box_red);

% Change the coordinates of bboxes to be coherent with
% other results from other sources (sel_search, etc.)
candidates.bboxes = [red_bboxes(:,2) red_bboxes(:,1) red_bboxes(:,4) red_bboxes(:,3)];

% Get the labels of leave regions that form each candidates
candidates.superpixels = f_lp;
if ~isempty(f_ms)
    candidates.labels = cands2labels(cand_labels,f_ms);
else
    candidates.labels = {1};
end

seg_obj.num_segs.final_num_segs = size(cand_labels,1);
seg_obj.timings.total_seg_time = toc(start_time) - total_load_time;
seg_obj.preload = preload;

fprintf('Total time %.2f secs\n', seg_obj.timings.total_seg_time);

% A. Humayun - Fall 2014
% you can also use masks = savecands2masks( candidates ); if loading from
% saved file

% Transform the results to masks
if compute_masks
    if ~isempty(f_ms)
        candidates.masks = cands2masks(cand_labels, f_lp, f_ms);
    else
        candidates.masks = true(size(f_lp));
    end
end

% storing the parameters
seg_obj.segm_params = struct;
seg_obj.segm_params.mode = mode;
seg_obj.segm_params.J_th = J_th;
seg_obj.segm_params.J_th_mex_box_red = J_th_mex_box_red;
seg_obj.segm_params.theta = theta;
seg_obj.segm_params.scales = scales;
seg_obj.segm_params.pareto_n_cands_fp = pareto_n_cands_fp;
seg_obj.segm_params.rf_regressor_fp = rf_regressor_fp;