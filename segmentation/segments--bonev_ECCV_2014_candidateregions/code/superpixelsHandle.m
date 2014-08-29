% Canidate Regions - version 1.0.  Code distributed under LGPL. 
% August 2014, University of California, Los Angeles
%
% Please, cite this work if you use the code: 
%   B. Bonev, A.L. Yuille, "A Fast and Simple Algorithm for Producing Candidate Regions", 
%   ECCV 2014, Zürich, Switzerland, 6th-12th September, 2014

classdef superpixelsHandle < handle
    properties
        ih;         % The image handler, not really needed here.
        RC;         % Rows, Cols
        segments;   % Rows x Cols labeled image starting from 1
        NR;         % Number of regions (of labels in "segments")
        st;         % regionprops statistics structure
        edgepriorW; % a weight, not used here
        AR;         % region areas in 1xNR 
        
        % Graph-related properties 
        A;          % Adjacency
        S;          % Edgeness sum between direct neighbours by imRAGb.m
        L;          % Length of borders between direct neighbours by imRAGb
        % Appearance properties
        EX;         % E[X];
        EX2;        % E[X^2]; also take into account the areas AR

        BB;         % 4 bounding box
        CC;         % 2 centroid
        
        s1;         % handler to the first level of segmentsHandle class
        s1corr;     % 1xNR containing the correspondences with first level
                    % It would usually be s1corr == 1:NR
        
        parts;      % helper for evaluation purposes
    end
    methods
        function obj = superpixelsHandle()
            obj.NR = 0;
        end
        
        %% Segment functions
        function obj = slic(obj, img, slicSize, slicRegu)
            %run /home/boyan/prog/ucla/vlfeat-0.9.16/toolbox/vl_setup
            obj.segments = vl_slic(img, slicSize, slicRegu)+1;
            %The following could be removed but then
            %the structure obj.s1corr would be sometimes necessary
            n = 0;
            diff = false;
            for i=1:max(obj.segments(:))
                pixels = (obj.segments == i);
                if any(pixels(:))
                    n = n + 1;
                    if diff
                        obj.segments(pixels) = n;
                    end
                else
                    diff = true;
                end
                
            end
            
            obj.setSegments(obj.segments);
        end
                
        function obj = regionprops(obj)
            obj.st = regionprops(obj.segments,...
                    'Centroid','Area',...
                    'BoundingBox','PixelIdxList');
            obj.AR = cat(2,obj.st(:).Area);
        end
        function obj = setSegments(obj, segments)
            obj.segments = uint32(segments);
            obj.NR = max(max(segments));
            obj.RC = size(segments);
        end
        
        %used in groundtruth segmentation handling
        function [newsegments correspondences] = ...
                 renumberSegments(obj, segments, objects)
            if ~exist('objects','var')
                objects = zeros(size(segments));
            end
            %assuming no more than 100 class labels in segments
            maxclasses = 100;
            newsegments = single(int32(segments)+int32(objects)*maxclasses);
            newnumber = 1;
            nr = max(max(newsegments));
            for r=1:nr
                found = find(newsegments==r);
                if isempty(found)
                    continue
                end
                correspondences(newnumber) = mod(r,maxclasses);
                newsegments(found) ...
                    = -newnumber; %negative label to distinguish it
                newnumber = newnumber + 1;
            end
            newsegments = -newsegments; %they had all negative labels
        end
        
        function initParts(obj) % For evaluation purposes, mainly.
            obj.parts = false([obj.NR, size(obj.segments)]);
            for s = 1:obj.NR
                part = false(size(obj.segments));
                part(obj.st(s).PixelIdxList) = true;
                obj.parts(s,:,:) = part;
            end
        end
        
        %% Graph functions
  
        function obj = adjacency(obj, edgeprior)

            [aristas astats] = imRAGb(obj.segments, edgeprior);
            
            obj.A = false(obj.NR);
            obj.S = double(zeros(obj.NR));
            for i=1:size(aristas,1)
                obj.A(aristas(i,1),aristas(i,2))=true;
                obj.A(aristas(i,2),aristas(i,1))=true;
                obj.L(aristas(i,1),aristas(i,2))=astats(i,1);
                obj.L(aristas(i,2),aristas(i,1))=astats(i,1);
                obj.S(aristas(i,1),aristas(i,2))=astats(i,2);
                obj.S(aristas(i,2),aristas(i,1))=astats(i,2);
            end        
            
        end
        
        %% Appearance-based similarity functions

        %ii is 1xN samples; returns NxD
        function [E1 E2] = getVectors(obj, ih, ii)
            nr = size(ii,2);
            %data = zeros(nr, 14); 
            E1 = zeros(nr,7);
            E2 = zeros(nr,7);
            for i = 1:nr
                idxs = obj.st(ii(i)).PixelIdxList(:);
                % Nsamples x Ndimensions
                D = [ih.labl(idxs),...
                     ih.laba(idxs),...
                     ih.labb(idxs),...
                     ih.gx(idxs),...
                     ih.gy(idxs),...
                     ih.gxx(idxs),...
                     ih.gyy(idxs),...
                     ];
                ni = 1/size(D,1);
                E1(i,:) = sum(D,1)*ni;
                E2(i,:) = sum(D.^2,1)*ni;
                %data(i,:) = [mean(D,1), std(D,1)];
            end
        end

%         %ii is 1xN samples; returns NxD
%         %ii should contain [i ineighbors]
%         function V = getNormalizedVectors(obj, ih, ii)
%             D = obj.getVectors(ih, 1:obj.NR);
%             means = mean(D,1);
%             stds  = std( D,1);
%             stds(stds==0) = 1; % Anyway mean would be 0, then 0/1 = 0
%         
%             ni = size(ii,2);
%             V = (D(ii,:) - repmat(means,ni,1)) ./ repmat(stds,ni,1);
%         end

        function allvectors(obj, ih)
            obj.ih = ih;
            %obj.V = obj.getNormalizedVectors(ih, 1:obj.NR);
            [obj.EX obj.EX2]  = obj.getVectors(ih, 1:obj.NR);
            obj.BB = cat(1,obj.st(1:obj.NR).BoundingBox);
            obj.CC = cat(1,obj.st(1:obj.NR).Centroid);
%            [obj.bEX obj.bEX2] = obj.getVectors(ih, 1:obj.NB, true);
        end
    end
end
