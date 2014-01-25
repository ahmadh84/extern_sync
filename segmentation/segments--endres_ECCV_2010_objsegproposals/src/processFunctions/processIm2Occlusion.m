function [bndinfo, pbim, gconf, bndinfo_all, pb_time] =  processIm2Occlusion(im, varargin)

%% Set parameters

%% Read image

if max(size(im))>640
  fprintf('Warning, this image is pretty big...\n');
  %im = imresize(im, 640/max(size(im)), 'bilinear');
end

%% Get occlusion info
[bndinfo, pbim, gconf, bndinfo_all, pb_time] = im2boundariesTopLevel(im);
gconf = single(gconf);
pbim = single(pbim);

%save(outname, 'bndinfo', 'pbim', 'gconf', 'bndinfo_all');
