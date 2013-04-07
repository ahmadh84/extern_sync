curr_dir = fileparts(which(mfilename));
utils_dir = fullfile(curr_dir, '..', '..');

addpath([curr_dir '/detector'])
addpath([curr_dir '/detector/util'])
addpath([curr_dir '/descriptor'])
addpath([curr_dir '/display'])

addpath(fullfile(utils_dir, 'graph_algos', 'matlab_bgl'))
addpath(fullfile(utils_dir, 'metrics', 'hua_pwmetric'))
addpath(fullfile(utils_dir, 'toolboxes', 'vlfeat', 'toolbox'))
addpath(fullfile(utils_dir, 'segmentation', 'boundaries+segments--arbelaez_PAMI_2010_bsr', 'lib'))
addpath(fullfile(utils_dir, 'graph_algos', 'spanning_tree'))
vl_setup
disp('BPLR ready.')