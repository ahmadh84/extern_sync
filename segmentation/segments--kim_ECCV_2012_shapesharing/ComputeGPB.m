function [ucm, ucm2, gpb, gpb_orient, t] = ComputeGPB(img_file)

disp('Compute gPb...');
t_pb = tic;
[gpb_orient, gpb, ~] = globalPb(img_file, '', 0.5);
t = toc(t_pb);

ucm2 = contours2ucm(gpb_orient, 'doubleSize');
ucm = ucm2(3:2:end, 3:2:end);