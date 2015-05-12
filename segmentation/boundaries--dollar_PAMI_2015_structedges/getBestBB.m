function [maxf, maxidx] = getBestBB( filepath, bbs,  segIds, edgeArr)
%BOTTLES Summary of this function goes here
%   Detailed explanation goes here

extern_path = fullfile(fileparts(which(mfilename)), '..', '..');
addpath(extern_path);
addpath(fullfile(extern_path, 'toolboxes/everingham_IJCV_2010_pascalvoc/VOCcode'));
addpath('/home/ahumayun/videovolumes/rigor_src/utils');

VOCinit;
[d f] = fileparts(filepath);

imshow('/data/images/everingham_IJCV_2010_pascalvoc/JPEGImages/2009_005302.jpg');

rec=PASreadrecord(sprintf(VOCopts.annopath,f));
anno = vertcat(rec.objects.bbox);

gtbbsz = anno;
gtbbsz(:,[3 4]) = gtbbsz(:,[3 4]) - gtbbsz(:,[1 2]) + 1;
oa = bbGt('compOas', bbs(:,1:4), gtbbsz);

[maxf,maxidx] = max(oa);

% save('bottles_ms1.mat', 'maxf','maxidx');
end

