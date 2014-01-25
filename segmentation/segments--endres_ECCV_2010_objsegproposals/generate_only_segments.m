function [proposal_data, superpixels, additional_info] = generate_only_segments(input)
%GENERATE_ONLY_SEGMENTS Summary of this function goes here
%   Detailed explanation goes here

% create multiple threads (set how many you have)
N_THREADS = 4;
if matlabpool('size') ~= N_THREADS
    matlabpool('open', N_THREADS);
end

function_root = which('generate_proposals.m');
function_root = fileparts(function_root);

proposals_startup;

col = load(fullfile(function_root, 'classifiers', 'colorClusters.mat'));
tex = load(fullfile(function_root, 'classifiers', 'textonClusters.mat'));
msclassifiers = load(fullfile(function_root, 'classifiers', 'msBgClassifiers.mat'));
load(fullfile(function_root, 'classifiers', 'bboxClassifier2.mat'), 'classifier_bbox');
load(fullfile(function_root, 'classifiers', 'subregionClassifier_mix.mat'), 'classifier');

additional_info.timing = struct;

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
[occ.bndinfo, occ.pbim, image_data.gconf, occ.bndinfo_all, pb_time] = ...
   processIm2Occlusion(im);                                        % random
[occ.pb1, occ.pb2, occ.theta] = getOrientedOcclusionProbs(occ.bndinfo_all);
bmaps = getOcclusionMaps(occ.bndinfo_all); 
occ.bmap = mean(bmaps,3); 


image_data.occ = occ;

additional_info.timing.pb = pb_time;
additional_info.timing.occlbndry_geomcont = toc(start_ob) - additional_info.timing.pb;
fprintf('Done (%f)\n', toc(start_ob));

%%%% Color + Texture codewords
start_ct = tic;
fprintf('------Quantize color/texture------\n');

[image_data.textonim, image_data.colorim] = processIm2ColorTexture(im, col, tex);
additional_info.timing.colortext = toc(start_ct);
fprintf('Done (%f)\n', toc(start_ct));

%%%% Probability of BG
fprintf('------Probability of BG------\n');
start_bg = tic;
[image_data.bg] = processIm2MsObjects(im, msclassifiers);          % random
additional_info.timing.prob_bg = toc(start_bg);
fprintf('Done (%f)\n', toc(start_bg));
fprintf('\nTotal time: %f\n', toc(start_image));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Segment level features %%%%%%%%%%%
fprintf('\n***Extracting Segment Features******\n');
start_seg = tic;
[region_data] = processData2RegionFeatures(image_data, classifier_bbox, classifier);
additional_info.timing.region_features = toc(start_seg);
fprintf('Done (%f)\n', toc(start_seg));

fprintf('\n***Proposing Regions******\n');
start_prop = tic;
[proposals proposal_data timings_prop] = proposeRegions(image_data, region_data);
fs = fields(timings_prop);
for idx = 1:length(fs);
    additional_info.timing.(fs{idx}) = timings_prop.(fs{idx});
end
fprintf('Done (%f)\n', toc(start_prop))

additional_info.timing.t_all = toc(start);
additional_info.num_seeds = length(region_data.regions);

superpixels = image_data.occ.bndinfo_all{1}.wseg;

end

