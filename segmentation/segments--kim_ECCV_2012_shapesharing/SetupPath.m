%% external codes
extern_src = '/home/ahumayun/videovolumes/extern_src/';

bplr_path = fullfile(extern_src, 'segmentation', 'segments--kim_CVPR_2011_bplr');
run(fullfile(bplr_path, 'bplr_startup.m'));

gc_path = fullfile(extern_src, 'optimization', 'maxflow--boykov_PAMI_2001_GCmex1.5');
addpath(gc_path)

gpb_path = fullfile(extern_src, 'segmentation', 'boundaries+segments--arbelaez_PAMI_2010_bsr');
addpath(genpath(gpb_path));

flann_path = fullfile(extern_src, 'machine_learning', 'flann', 'build', 'src', 'matlab');
%addpath(flann_path);
flann_path = fullfile(extern_src, 'machine_learning', 'flann', 'src', 'matlab');
addpath(flann_path);

curr_path = fileparts(which(mfilename));
external_code_path = fullfile(curr_path, 'external_code');

flann_path = fullfile(external_code_path, 'flann-1.21-src-64', 'build', 'matlab');
addpath(flann_path);

addpath(fullfile(curr_path, 'external_code'));

%% mex files
addpath(fullfile(curr_path, 'mex'));

%% pascal
addpath(fullfile(curr_path, 'VOCcode'));
