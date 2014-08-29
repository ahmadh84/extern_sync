% Canidate Regions - version 1.0.  Code distributed under LGPL. 
% August 2014, University of California, Los Angeles
%
% Please, cite this work if you use the code: 
%   B. Bonev, A.L. Yuille, "A Fast and Simple Algorithm for Producing Candidate Regions", 
%   ECCV 2014, Zürich, Switzerland, 6th-12th September, 2014
function [segh] = hdsegments3(image, doplot, model, edgesWeight, slicSize, slicRegu, ratioMerge, nfewsegments)
% Input:    image in RGB or image path
% Optional: doplot is true by default to show plots during segmentation
%           edgeprior grayscale image with edges
%           edgepriorWeight from 0 to approximately 10. default is 1
% Output:   segh is an array of segmentsHandle handles to objects
%           segh(i).segments contains the segments map of level i
% Boyan Bonev - bonev@ucla.edu - April 2013

    
    if isempty(edgesWeight)
        edgesWeight = 1.3;          
    end
    if isempty(model),  
        ih = imageHandle(image); 
    else
        %%%edgeprior = imdilate(single(edgeprior),strel('square',2))/255; 
        %image2 = uint8(255*single_correction(double(image)/255, 0.4));
        ih = imageHandle(image, edgesDetect(image,model));   %Piotr Dollar          
    end

    %First level segmentation with SLIC
    %slicSize = floor(prod(ih.RC)*0.0001); % 30 pixels for 640x480
    %if slicSize<10, slicSize=10; end

    %slicRegu = 0.1;
    s0 = superpixelsHandle();
    s0.slic(ih.sing, slicSize, slicRegu);
    s0.edgepriorW = edgesWeight;
    s0.regionprops();
    s0.adjacency(ih.edg); 
    s0.allvectors(ih);

    sh = segmentsHandle(1);
    sh.setFirstLevel(s0);
    segh(1) = sh;

    %    load('../Piotr/SketchTokens-master/models/forest/modelSmall.mat');
    %    st = stDetect( image, model );
    %    edgeprior = stToEdges( st, 1 ) * 255;
    
    if doplot, 
        figure, imshow(uint8(s0.edgeprior*255));
        figure, imshow(ih.img); 
    end
    
    for level = 1:100 % A level in the segmentation hierarchy
        [sh2] = pairwiseMerge(sh, ratioMerge, nfewsegments, doplot);

        segh(level+1) = sh2; %Note that they are handlers (pointers)
        sh = sh2;
        
        if sh2.NR <= 1 || sh2.NR == segh(level).NR 
            break
        end
    end
    segh(end).parts = 1:segh(end).NR;

end

%Mergepercentage should be between 0 and 1, better less than 0.5
%Highlevel indicates at which level in the hierarchy k should becom 1,
%because segments are already too big and makes no sense merging to a 
%second nn.
function [sh2] = pairwiseMerge(sh, mergepercentage, fewsegments, doplot)
    
    %Compute distances
    W = zeros(sh.NR);  % distances graph
    RW = zeros(sh.NR); % 1/distances graph for page ranking
    NG = zeros(sh.NR); % on each row, the region's neighbors
    
    if sh.NR<fewsegments 
        k = 1;
    else
        k = 2;
    end
    sh.precalcW(k); %precalculates all distances that will be necessary
    for r = 1:sh.NR
%        [dists neighbors bestv] = sh.dist2neighbors(r, k); %first the good matches
        
        ineighbors = find(sh.W(r,:)~=-1);
        dists = sh.W(r, ineighbors);
        if sh.NR <= fewsegments %obj.level > 27
            dists = dists.*(0.05 + 0.95 ./ (sh.L(r,ineighbors).*max(sh.L(:))) .* (sh.AR(ineighbors)./sh.AR(r)).^2 );
        end
        [~,ii] = sort(dists,2,'ascend'); %first the good matches
        dists = dists(ii);
        neighbors = ineighbors(ii);
        
        %dists = dists + 0.00000000001*rand(size(dists)); %just in case of regularities
        dists(dists==0) = 0.00000000001; % Residual for not dividing by 0
        
        NG(r,1:size(neighbors,2)) = neighbors;
        
        W(r,neighbors) = dists; %dissimilarity
        W(r,r) = 0;
        
        %Modify weights for ranking here. However manipulation is not good.
        RW(r,neighbors) = 1./dists; %links from r --> to nns
        RW(r,r) = 0;
        
    end
    
    RW = RW.*sh.L; %weight similarity by length of common border

    % Rank nodes
    matriz = (RW ./repmat(sum(RW,2),1,sh.NR))';
    %matriz = (matriz + (matriz)')/2;
    p=rank( matriz ,0.85, 0.0001);
    %p = rand(sh.NR,1);
    
    sh.rank = [0.7, 0.3]*[sh.rank; p']; %weighted mean
    [~, sorted] = sort(sh.rank,2,'descend'); %higher ranks first
    
    % Merge first nodes in ranking (most 'desired' ones, or most 'popular').
    associated = [];
    merge1=zeros(sh.NR,1);
    merge2=zeros(sh.NR,1);
    i = 1;
    for r = sorted
        closer = NG(r,1);
                
        if ~isempty(find(associated==r, 1)) || ~isempty(find(associated==closer, 1)) % || closer==0
            continue
        end
        
        %filter bad neighbors
        closernn = NG(closer,find(NG(closer,:)));
        goodnn = closernn(1:ceil(end/2));
        if isempty(find(goodnn==r)) || W(closer,r)>W(r,closer)*1.15
            %The first condition means that I am not among the
            % first options of the neighbor that I've chosen as 1st opt.
            %The second condition means that the distance from
            % my neighbor to me can only be slightly larger than 
            % the distance from me to the neighbor. I don't want him
            % to dist much more than I do from him, but as I (r) am
            % choosing, no problem if I dists much more than he does:
            % he likes me much, I don't like it that much, but I
            % choose him. This is necessary, otherwise it gets stuck.
            continue;
        end
        
        associated = [associated closer r]; 
        %Or, more conservatively, also disable touching my neighbours:
        %associated = [associated  r  sh.getKstepsNeighbours(r,1) closer sh.getKstepsNeighbours(closer,1)]; aprop=0.6;
        
        merge1(i) = r; merge2(i) = closer; %and merge later
        i=i+1;
        
        if size(associated,2) >= round(sh.NR*mergepercentage)
           break;
        end
        
    end
    sh.associated = associated;
    sh2 = segmentsHandle(sh.level + 1);
    sh2.setNextLevel(sh);
    sh2.mergeSegments(merge1(1:i-1),merge2(1:i-1)); %Merging finished.
    sh.parts = [merge1(1:i-1); merge2(1:i-1)]';
    
    if doplot
        figure, subplot(2,2,1);
        title(['Solidity level ' int2str(sh.level)]);
        subplot(2,2,2);
        imagesc(sh.getSegmentsPlot()), title(['Segmentation level ' int2str(sh.level)]);
        subplot(2,2,3);
        rankplot = sh.plotGraphsHeat(sh.rank);
        imagesc(rankplot), title(['Rank ' int2str(sh.level)]);
        subplot(2,2,4);
        imagesc(sh.maskRegions(sh.parts)), title('Output segments');
    end
    
end

