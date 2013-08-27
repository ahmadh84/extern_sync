curr_dir = fileparts(which(mfilename));

addpath(genpath(fullfile(curr_dir, 'src')));
addpath(genpath(fullfile(curr_dir , 'external')));
addpath(fullfile(curr_dir, '..', 'segbench', 'lib', 'matlab'));
addpath(fullfile(curr_dir , '..', 'segments--felzenszwalb_IJCV_2004_superpixels'));
addpath(genpath(fullfile(curr_dir , '..', '..', 'toolboxes', 'lightspeed')));
addpath(genpath(fullfile(curr_dir , '..', '..', 'metrics', 'hua_pwmetric')));
addpath(fullfile(curr_dir , '..', '..', 'image_processing', 'bresenham_mex'));
addpath(fullfile(curr_dir , '..', '..', 'optimization', 'maxflow--boykov_PAMI_2004_maxflow', 'mex_maxflow'));
