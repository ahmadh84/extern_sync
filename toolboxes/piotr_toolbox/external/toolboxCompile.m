% Compiles all mex routines that are part of toolbox.
%
% USAGE
%  toolboxCompile
%
% INPUTS
%
% OUTPUTS
%
% EXAMPLE
%
% See also
%
% Piotr's Image&Video Toolbox      Version 3.22
% Copyright 2013 Piotr Dollar.  [pdollar-at-caltech.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Simplified BSD License [see external/bsd.txt]

% compile options including openmp support for C++ files
opts = {'-output'};
if(exist('OCTAVE_VERSION','builtin')), opts={'-o'}; end
if( ispc ), optsOmp={'OPTIMFLAGS="$OPTIMFLAGS','/openmp"'}; else
  optsOmp={'CXXFLAGS="\$CXXFLAGS','-fopenmp"'};
  optsOmp=[optsOmp,'LDFLAGS="\$LDFLAGS','-fopenmp"'];
end

% list of files (missing /private/ part of directory)
fs={'channels/convConst.cpp', 'channels/gradientMex.cpp',...
  'channels/imPadMex.cpp', 'channels/imResampleMex.cpp',...
  'channels/rgbConvertMex.cpp', 'classify/binaryTreeTrain1.cpp', ...
  'classify/fernsInds1.c', 'classify/forestFindThr.cpp',...
  'classify/forestInds.cpp', 'classify/meanShift1.c',...
  'detector/acfDetect1.cpp', 'images/assignToBins1.c',...
  'images/histc2c.c', 'images/imtransform2_c.c', ...
  'images/nlfiltersep_max.c', 'images/nlfiltersep_sum.c', ...
  'videos/ktComputeW_c.c', 'videos/ktHistcRgb_c.c', ...
  'videos/opticalFlowHsMex.cpp' };
n=length(fs); useOmp=zeros(1,n); useOmp([6 9])=1;

% compile every funciton in turn (special case for dijkstra)
disp('Compiling Piotr''s Toolbox.......................');
rd=fileparts(mfilename('fullpath')); rd=rd(1:end-9); tic;
try
  for i=1:n
    [d,f1,e]=fileparts(fs{i}); f=[rd '/' d '/private/' f1];
    if(useOmp(i)), optsi=[optsOmp opts]; else optsi=opts; end
    fprintf(' -> %s\n',[f e]); mex([f e],optsi{:},[f '.' mexext]);
  end
  d=[rd '/matlab/private/']; fprintf(' -> %s\n',[d 'dijkstra1.cpp']);
  mex([d 'fibheap.cpp'], [d 'dijkstra1.cpp'], '-largeArrayDims', ...
    opts{:}, [d 'dijkstra1.' mexext]);
catch ME
  fprintf(['C++ mex failed, likely due to lack of a C++ compiler.\n' ...
    'Run ''mex -setup'' to specify a C++ compiler if available.\n'...
    'Or, one can specify a specific C++ explicitly (see mex help).\n']);
end
disp('..................................Done Compiling'); toc;
