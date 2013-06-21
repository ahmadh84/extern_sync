function dummy = example(IFolder,imgFile1);

%% Compute global Pb and hierarchical segmentation for an example image.
addpath(fullfile(pwd,'lib'));

%% 1. compute globalPb
imgFilePath = IFolder;
imgFileEnd = '.ppm';
imgFile = sprintf('%s%s%s',imgFilePath,imgFile1,imgFileEnd);
outFile = 'data/101087_gPb.mat';
gPb_orient = globalPb(imgFile, outFile);

%% 2. compute Hierarchical Regions

% for boundaries
ucm = contours2ucm(gPb_orient, 'imageSize');

% for regions
ucm2 = contours2ucm(gPb_orient, 'doubleSize');

%% 3. usage example
% convert ucm to the size of the original image
ucm = ucm2(3:2:end, 3:2:end);

% write probability boundary set
 imwrite(ucm,sprintf('%s%s_bdrySET%s',imgFilePath,imgFile1,imgFileEnd));

% ---------------------------------------
% superpixels to piecewise constant image
% ---------------------------------------
k = 0.05;
bdry = (ucm >= k);
labels2 = bwlabel(ucm2 <= k);
labels = labels2(2:2:end, 2:2:end);
Sp2 = labels;
m = max(max(Sp2));
d = floor(32768/m);
A(:,:,1) = Sp2*d;
A(:,:,2) = Sp2*d;
A(:,:,3) = Sp2*d;
A(:,:,1) = floor(A(:,:,1) / (32*32));
A(:,:,2) = mod(A(:,:,2), 32*32);
A(:,:,2) = floor(A(:,:,2) / 32);
A(:,:,3) = mod(A(:,:,3), 32);
A = A * 8;
%imwrite(bdry,sprintf('%s%s_bdry0%s',imgFilePath,imgFile1,imgFileEnd));
imwrite(uint8(A),sprintf('%s%s_seg0%s',imgFilePath,imgFile1,imgFileEnd));

k = 0.1;
bdry = (ucm >= k);
labels2 = bwlabel(ucm2 <= k);
labels = labels2(2:2:end, 2:2:end);
Sp2 = labels;
m = max(max(Sp2));
d = floor(32768/m);
A(:,:,1) = Sp2*d;
A(:,:,2) = Sp2*d;
A(:,:,3) = Sp2*d;
A(:,:,1) = floor(A(:,:,1) / (32*32));
A(:,:,2) = mod(A(:,:,2), 32*32);
A(:,:,2) = floor(A(:,:,2) / 32);
A(:,:,3) = mod(A(:,:,3), 32);
A = A * 8;
%imwrite(bdry,sprintf('%s%s_bdry1%s',imgFilePath,imgFile1,imgFileEnd));
imwrite(uint8(A),sprintf('%s%s_seg1%s',imgFilePath,imgFile1,imgFileEnd));

k = 0.15;
bdry = (ucm >= k);
labels2 = bwlabel(ucm2 <= k);
labels = labels2(2:2:end, 2:2:end);
Sp2 = labels;
m = max(max(Sp2));
d = floor(32768/m);
A(:,:,1) = Sp2*d;
A(:,:,2) = Sp2*d;
A(:,:,3) = Sp2*d;
A(:,:,1) = floor(A(:,:,1) / (32*32));
A(:,:,2) = mod(A(:,:,2), 32*32);
A(:,:,2) = floor(A(:,:,2) / 32);
A(:,:,3) = mod(A(:,:,3), 32);
A = A * 8;
%imwrite(bdry,sprintf('%s%s_bdry2%s',imgFilePath,imgFile1,imgFileEnd));
imwrite(uint8(A),sprintf('%s%s_seg2%s',imgFilePath,imgFile1,imgFileEnd));

