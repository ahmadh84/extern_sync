function varargout = imRAGb(img, varargin)
%IMRAG Region adjacency graph of a labeled image
% by David Legland. See BSD License. 

%Copyright (c) 2012, INRA
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are
%met:
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the distribution
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%POSSIBILITY OF SUCH DAMAGE.

%   Modification by Boyan Bonev, bonev@ucla.edu, March 2013:
%   Apart from edges returns boundary size and two additional matrices
%   with labeled boundaries for each region: from and to:
%   [edges edgestats bound_from bound_to] = imRAGb(img, appearanceedges);
%   edges is Nx4, bound_xx are of img size.
%   edges has:        edgestats has:
%   i  j              edge_size mean_value
%   ...
%   where mean_value are the mean values of "appearanceedges" under the
%   boundaries of the neighbors.
%     Example with regions 1,2 and 4:
%     a =
%          1     1     2     2
%          1     1     2     2
%          2     2     4     4
%          2     2     4     4
%     e = % some intensities whose means matter on the edges
%                0           0           0           0
%                0           6          50        1000
%                0           0           0           0
%             1000           0           0           0
% 
%     [aristas, stats, bf, bt] = imRAGb(a,e)
%     aristas = % see rows:
%          1     2   
%          2     4   
%     stats =   % see rows:
%          4     8
%          4   150
%     bf =
%          0     1     2     0
%          1     1     2     2
%          2     2     4     4
%          0     2     4     0
%     bt =
%          0     2     1     0
%          2     2     4     4
%          1     4     2     2
%          0     4     2     0
% Note that edges are double: one pixel on each region

% size of image
dim = size(img);

% output matrices with boundaries
% bnd_from = single(zeros(dim));
% bnd_to   = single(zeros(dim));

if ~isempty(varargin)
    onboundaryedgeness = varargin{1};
else
    onboundaryedgeness = zeros(dim);
end


%% First direction of 2D image

% identify transitions
[i1 i2] = find(img(1:end-1,:) ~= img(2:end, :));

% remove values close to border
inds = i1 >= dim(1);
i1(inds) = [];
i2(inds) = [];

% get values of consecutive changes
plist1 = sub2ind(dim, i1, i2);
plist2 = sub2ind(dim, i1+1, i2);
val1 = img(plist1);
val2 = img(plist2);

% find all changes not involving background
inds = val1 ~= 0 & val2 ~= 0 & val1 ~= val2;

edges = [val1(inds) val2(inds)];
plist = [plist1(inds) plist2(inds)];

%% Second direction of 2D image

% identify transitions
[i1 i2] = find(img(:, 1:end-1) ~= img(:, 2:end));

% remove values close to border
inds = i2 >= dim(2);
i1(inds) = [];
i2(inds) = [];

% get values of consecutive changes
plist1 = sub2ind(dim, i1, i2);
plist2 = sub2ind(dim, i1, i2+1);
val1 = img(plist1);
val2 = img(plist2);

% find all changes not involving background
inds = val1 ~= 0 & val2 ~= 0 & val1 ~= val2;

edges = [edges; [val1(inds) val2(inds)]];
plist = [plist; [plist1(inds) plist2(inds)]];

%% Summarize the results
uedges= unique(sortrows(sort(edges, 2)),'rows');
summs = zeros(size(uedges,1),1);
lengs = double(zeros(size(uedges,1),1));
for i=1:size(uedges,1)
    %vertical changes:
    found1 = find( edges(:,1)==uedges(i,1) & edges(:,2)==uedges(i,2) );
    %horizontal changes:
    found2 = find( edges(:,1)==uedges(i,2) & edges(:,2)==uedges(i,1) );
    %found1 and found2 refer to different changes (edges) but may describe
    %repeated pixels in the image
    uind = [sets.unique([plist(found1,1); plist(found1,2)]);  ...
            sets.setdiff([plist(found2,1); plist(found2,2)], sets.intersect(plist(found1,2),plist(found2,1))  )   ];
    lengs(i) = numel(uind);
    summs(i) = sum(onboundaryedgeness(uind));
end

edges = uedges;


%% Output processing

if nargout >= 1
    varargout{1} = edges;
end
if nargout >=2
    varargout{2} = [double(lengs)*0.5 summs*0.5];
end

end


