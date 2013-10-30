function [proposal_data, superpixels, total_time] = generate_only_segments(input)
%GENERATE_ONLY_SEGMENTS Summary of this function goes here
%   Detailed explanation goes here

function_root = which('generate_proposals.m');
function_root = fileparts(function_root);

proposals_startup;

col = load(fullfile(function_root, 'classifiers', 'colorClusters.mat'));
tex = load(fullfile(function_root, 'classifiers', 'textonClusters.mat'));
msclassifiers = load(fullfile(function_root, 'classifiers', 'msBgClassifiers.mat'));
load(fullfile(function_root, 'classifiers', 'bboxClassifier2.mat'), 'classifier_bbox');
load(fullfile(function_root, 'classifiers', 'subregionClassifier_mix.mat'), 'classifier');


if(isstr(input))
   im = im2double(imread(input));
else
   im = im2double(input); % Just to make input consistent
end

start = tic;
start_image = start;
fprintf('***Extracting image level features******\n');

%%%% Occlusion and Geometric context
fprintf('------Occlusion boundaries + Geometric context------\n')
start_ob = tic;
[occ.bndinfo, occ.pbim, image_data.gconf, occ.bndinfo_all] = ...
   processIm2Occlusion(im);                                        % random
[occ.pb1, occ.pb2, occ.theta] = getOrientedOcclusionProbs(occ.bndinfo_all);
bmaps = getOcclusionMaps(occ.bndinfo_all); 
occ.bmap = mean(bmaps,3); 


image_data.occ = occ;

fprintf('Done (%f)\n', toc(start_ob));

%%%% Color + Texture codewords
start_ct = tic;
fprintf('------Quantize color/texture------\n');

[image_data.textonim, image_data.colorim] = processIm2ColorTexture(im, col, tex);
fprintf('Done (%f)\n', toc(start_ct));

%%%% Probability of BG
fprintf('------Probability of BG------\n');
start_bg = tic;
[image_data.bg] = processIm2MsObjects(im, msclassifiers);          % random
fprintf('Done (%f)\n', toc(start_bg));
fprintf('\nTotal time: %f\n', toc(start_image));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Segment level features %%%%%%%%%%%
fprintf('\n***Extracting Segment Features******\n');
start_seg = tic;
[region_data] = processData2RegionFeatures(image_data, classifier_bbox, classifier);
fprintf('Done (%f)\n', toc(start_seg));

fprintf('\n***Proposing Regions******\n');
start_prop = tic;
[proposals proposal_data] = proposeRegions(image_data, region_data);

fprintf('Done (%f)\n', toc(start_prop))

total_time = toc(start);

superpixels = image_data.occ.bndinfo_all{1}.wseg;

end

