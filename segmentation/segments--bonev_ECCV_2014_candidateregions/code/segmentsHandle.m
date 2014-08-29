% Canidate Regions - version 1.0.  Code distributed under LGPL. 
% August 2014, University of California, Los Angeles
%
% Please, cite this work if you use the code: 
%   B. Bonev, A.L. Yuille, "A Fast and Simple Algorithm for Producing Candidate Regions", 
%   ECCV 2014, Zürich, Switzerland, 6th-12th September, 2014

classdef segmentsHandle < handle
    properties
        s0;         % Level0 superpixels (handle)
        shprev;     % Previous level 
        level;      % Integer annotation of the segmentation level
        NR;         % Number of regions (of labels in "segments")
        edgepriorW; % Weight for edge prior, from 0 to 1
        AR;         % region areas in 1xNR 
        sqrtAR;     % sqrt(AR)
        EX;         % E[X];
        EX2;        % E[X^2]; also take into account the areas AR
        depths;     % Weights used when taking means of regions
        rank;       % page rank history means, for stabilization
        
        BB;         % NRx4 bounding box
        CC;         % NRx2 centroid
        
        % Normalization variables (for the vector: EX, STD, BB, CC)
        means;
        stds;
        
        % Graph-related properties 
        A;          % Adjacency
        S;          % Edgeness sum between direct neighbours by imRAGb.m
        L;          % Length of borders between direct neighbours by imRAGb
        E;          % Edgeness between direct neighbors:
                    % intended use: E = S./L;
        W;          % Weights matrix: the bidirectional distances
                    % Between first and second neighbors, using edgeness.
                    % to exclude edgeness: W2=obj.W-obj.edgepriorW*(~obj.A + obj.E)
        
        % Structure relating segments with superpixels in s0
        hass;       % NR x s0.NR with arrays on rows, indicating s0 indexes
        hasn;       % NR x 1 : how many non-zeros has each row in hass
        % Structure relating segments segments with previous level segmts
        hads;       % NR x 2 segment indexes from previous level
                    % can have between 1 or 2 segments associated.
        
        % Higher level statistics (Optional)
        parts;      % Those object-like segments which will be merged next
        associated; % All segments (including parts) which will be merged
        complexity; % The complexity of the graph
        
        % For segments selection
        heterogen;  % 1xNR temporary structure, -1 if a cell has no value
    end
    methods
        %% Construction-related functions
        function obj = segmentsHandle(level)
            obj.level = level;
            obj.NR = 0;
            obj.complexity = 0;
            obj.W = [];
            obj.heterogen = [];
        end
        
        function setFirstLevel(obj, s0)
            obj.s0   = s0; obj.s0.s1 = obj;
            obj.shprev=s0;
            idx=find(cat(1,s0.st(:).Area)>0); %nonzero superpixels
            obj.NR   = size(idx,1);
            obj.s0.s1corr  = zeros(1,s0.NR); obj.s0.s1corr(idx) = 1:obj.NR;
            obj.edgepriorW = s0.edgepriorW;
            obj.hass = zeros(obj.NR,s0.NR);
            obj.hass(:,1) = idx;
            obj.hasn = ones(obj.NR,1);
            obj.hads = zeros(obj.NR,2);
            obj.hads(:,1)   = [1:obj.NR]';
            obj.A    = s0.A(idx,idx);
            obj.S    = s0.S(idx,idx);
            obj.L    = s0.L(idx,idx);
            nonzeros = obj.S~=0;
            obj.E    = zeros(obj.NR);
            obj.E(nonzeros) = obj.S(nonzeros)./obj.L(nonzeros);
            obj.AR   = cat(2,obj.s0.st(idx).Area);
            obj.sqrtAR=sqrt(obj.AR);
            obj.EX   = s0.EX(idx,:);
            obj.EX2  = s0.EX2(idx,:);
            obj.depths = ones(1,obj.NR);
            obj.BB   = s0.BB(idx,:);
            obj.CC   = s0.CC(idx,:);
            obj.rank = zeros(1,obj.NR);
            obj.initVectorNormalization();
            obj.W    = -1*ones(obj.NR); obj.W(logical(eye(obj.NR))) = 0;
        end
        
        function setNextLevel(obj, sh)
            obj.s0   = sh.s0;
            obj.shprev=sh;
            obj.NR   = sh.NR;
            obj.edgepriorW = sh.edgepriorW;
            obj.hass = sh.hass;
            obj.hasn = sh.hasn;
            obj.hads = zeros(obj.NR,2);
            obj.hads(:,1)   = [1:obj.NR]';
            obj.A    = sh.A;
            obj.S    = sh.S;
            obj.L    = sh.L;
            obj.E    = sh.E;
            obj.AR   = sh.AR;
            obj.sqrtAR=sqrt(obj.AR);
            obj.EX   = sh.EX;
            obj.EX2  = sh.EX2;
            obj.depths=sh.depths;
            obj.BB   = sh.BB;
            obj.CC   = sh.CC;
            obj.rank = sh.rank;
            obj.means= sh.means;
            obj.stds = sh.stds;
        end
        
        %% Graph functions
        function nn = getKstepsNeighbours(obj, i,Kstepsneighbors) % 1 or 2
            nn = find(obj.A(i,:));
            if Kstepsneighbors==2
                idx = any(obj.A(nn,:),1);
                idx(nn) = true; idx(i) = false;
                nn = find(idx);
            end
            %nn = setdiff(nn,i);
            %idx = false(1,obj.NR);
            %idx(nn) = 1; idx(i) = 0;
            %nn = find(idx);
            
        end
        
        %maxdepth has to be < obj.level
        %nodes, levels and updates are row vectors
        %updated is logical. use: nodes(updated) levels(updated)
        function [nodes levels updated] = getBranch(obj, node, maxdepth)
            nodes = [node];
            levels = [obj.level];
            updated = [1];
            nlastlevel = 1;
            sh = obj;
            for d = obj.level:-1:obj.level-maxdepth+1
                if sh.level <= 1
                    break
                end
                shhads = sh.hads(nodes(end-nlastlevel+1:end),:);
                newnodes = shhads(:)'; 
                newnodes = newnodes(newnodes~=0);
                a = logical(shhads(:,2))'; 
                newupdates = [a a(a)];
                nodes = [nodes, newnodes];
                updated = [updated, newupdates];
                nlastlevel = size(newnodes,2);
                levels = [levels, (sh.level-1)*ones(1,nlastlevel)];
                sh = sh.shprev;
            end
        end
        
        function [nns nnlevels] = getMultilevelNN(obj, node, maxdepth, mincontact)
            firstlevelnns = find(obj.A(node,:)); 
            firstlevelnns = firstlevelnns(...
                obj.L(node,firstlevelnns)./obj.sqrtAR(firstlevelnns) ...
                >= mincontact &...
                obj.L(node,firstlevelnns)/obj.sqrtAR(node) ...
                >= mincontact );

            nns = firstlevelnns;
            nnlevels = [obj.level*ones(size(firstlevelnns))];
            ispart = true(size(firstlevelnns));
            nodes = obj.hads(node,:); nodes = nodes(nodes~=0);
            sh = obj.shprev;
            for l = obj.level-1:-1:max(obj.level-maxdepth,1)
                A = sh.A;
                anyA = any(A(nodes,:),1);
                nextlevelnns = find(anyA);
                mask = false(1,sh.NR); mask(sh.parts) = true;
                mask(nodes) = false;
                nns = [nns, nextlevelnns];
                nnlevels = [nnlevels, l*ones(size(nextlevelnns))];
                ispart = [ispart, mask(nextlevelnns)];
                nodes = reshape(sh.hads(nodes,:),1,[]); nodes = nodes(nodes~=0);
                sh = sh.shprev;
            end
            
            nns = nns(ispart);
            nnlevels = nnlevels(ispart);
            
        end        
        
        function [nns nnlevels] = getMultilevelNNold(obj, node, maxdepth, mincontact)
            [n0 l0 u0] = obj.getBranch(node, maxdepth);
            firstlevelnns = find(obj.A(node,:)); 

            firstlevelnns = firstlevelnns(...
                obj.L(node,firstlevelnns)./obj.sqrtAR(firstlevelnns) ...
                >= mincontact &...
                obj.L(node,firstlevelnns)/obj.sqrtAR(node) ...
                >= mincontact );
            
            nns = [firstlevelnns];
            nnlevels = [obj.level*ones(size(firstlevelnns))];
            for nn0 = firstlevelnns;
                [ni li ui] = obj.getBranch(nn0, maxdepth);
                
                %at each level, check for common neighbors
                sh = obj.shprev;
                for level = obj.level-1:-1:obj.level-maxdepth
                    if level <= 1 
                        break
                    end
                    n0level = n0(l0==level);
                    %nilevel = ni(li==level & ui); %only updated ones
                    nilevel = ni(li==level);
                    nilevel = sets.intersect(nilevel,sh.parts)'; %only parts!
                    nn = [];
                    for nil = nilevel
                        
                        for n0l = n0level
                            nn  = [nn sets.intersect(sets.intersect(...
                                   [n0l find(sh.A(n0l,:))], ...
                                   [nil find(sh.A(nil,:))])',...
                                              nilevel)'...
                                  ];
                        end
                    end
                    nn = sets.unique(nn)';
                    nns = [nns nn];
                    nnlevels = [nnlevels level*ones(size(nn))];
                    sh = sh.shprev;
                end
            end
            %clean repeated
            unns = [];
            unnl = [];
            for level = obj.level:-1:obj.level-maxdepth
                un = sets.unique(nns(nnlevels==level))';
                unns = [unns un];
                unnl = [unnl level*ones(size(un))];
            end
            nns = unns;
            nnlevels = unnl;
            if isempty(nns)
                nns = [];
                nnlevels = [];
            end
        end
              
        % Indexes cannot be repeated in i and j. If repeated,
        % then the first merge will be valid but the next will remain as
        % a separate region.
        % Arrays ii and jj are the same size: nx1
        function mergeSegments(obj,ii,jj)
            for c=1:size(ii,1) %TODO: remove for by constraining ii,jj
                i=ii(c);
                j=jj(c);
                ni=obj.hasn(i);
                nj=obj.hasn(j);
                obj.hass(i,ni+1:ni+nj) = obj.hass(j,1:nj); %append j's to i row
                obj.hasn(i) = ni+nj;
                obj.hasn(j) = 0; %and later remove that row from hass and hasn
            end
            obj.hads(ii,2) = obj.hads(jj,1);
            De = obj.AR([ii,jj]).*obj.depths([ii,jj]); %two cols for i,j
            %De = obj.AR([ii,jj]); %No depths
            sumDe = sum(De,2);
            De = De ./ [sumDe, sumDe]; % norm to sum 1 on rows
            d = size(obj.EX,2);

            obj.EX(ii,:)  = ( (De(:,1)*ones(1,d)) .* obj.EX(ii,:) ) ...
                          + ( (De(:,2)*ones(1,d)) .* obj.EX(jj,:) );
            obj.EX2(ii,:) = (((De(:,1)*ones(1,d)) .* sqrt(obj.EX2(ii,:)))...
                          +  ((De(:,2)*ones(1,d)) .* sqrt(obj.EX2(jj,:)))...
                            ).^2;
                        
            We = obj.AR([ii,jj]); %two cols for i,j
            sumWe = sum(We,2);
            We = We ./ [sumWe, sumWe]; % norm to sum 1 on rows
            obj.depths(ii) = sum(We.*(obj.depths([ii,jj])+1),2);
            %Don't use obj.depths from this point on (in this function)
            obj.rank(ii)   = sum(We.* obj.rank([ii,jj]),2);
            obj.CC(ii,:)   = [... % weighted mean on x and on y
                  sum(We .* [obj.CC(ii,1), obj.CC(jj,1)], 2),... % CCx
                  sum(We .* [obj.CC(ii,2), obj.CC(jj,2)], 2)];   % CCy
            bbulx  = [obj.BB(ii,1) obj.BB(jj,1)];
            bbuly  = [obj.BB(ii,2) obj.BB(jj,2)];
            bbminx = min(bbulx,[],2); %col
            bbminy = min(bbuly,[],2); %col
            bbmaxx = max([bbulx(:,1)+obj.BB(ii,3), bbulx(:,2)+obj.BB(jj,3)],[],2);
            bbmaxy = max([bbuly(:,1)+obj.BB(ii,4), bbuly(:,2)+obj.BB(jj,4)],[],2);
            obj.BB(ii,:) = [bbminx, bbminy, bbmaxx-bbminx, bbmaxy-bbminy];
            obj.AR(ii) = sum(obj.AR([ii,jj]),2);
            %Don't read obj.AR from this point on (in this function)
            
            %Now, graph transformations:
            obj.A(ii,:)  = obj.A(ii,:)  | obj.A(jj,:);
            obj.A(ii,ii) = obj.A(ii,ii) | obj.A(ii,jj);
            obj.A(:,ii)  = obj.A(ii,:)';
            obj.A((ii-1)*(obj.NR+1)+1)  = 0;
            obj.S(ii,:)  = obj.S(ii,:)  + obj.S(jj,:);
            obj.S(ii,ii) = obj.S(ii,ii) + obj.S(ii,jj);
            obj.S(:,ii)  = obj.S(ii,:)';
            obj.S((ii-1)*(obj.NR+1)+1)  = 0;
            obj.L(ii,:)  = obj.L(ii,:)  + obj.L(jj,:);
            obj.L(ii,ii) = obj.L(ii,ii) + obj.L(ii,jj);
            obj.L(:,ii)  = obj.L(ii,:)';
            obj.L((ii-1)*(obj.NR+1)+1)  = 0;
            
            %remove zeroed rows:
            nonzero = obj.hasn>0;
            obj.hass = obj.hass(nonzero,:);
            obj.hasn = obj.hasn(nonzero);
            obj.hads = obj.hads(nonzero,:);
            obj.NR = size(obj.hasn,1);
            obj.CC = obj.CC(nonzero,:);
            obj.BB = obj.BB(nonzero,:);
            obj.EX = obj.EX(nonzero,:);
            obj.EX2 = obj.EX2(nonzero,:);
            obj.AR = obj.AR(nonzero);
            obj.sqrtAR = sqrt(obj.AR);
            obj.A = obj.A(nonzero,nonzero);
            obj.S = obj.S(nonzero,nonzero);
            obj.L = obj.L(nonzero,nonzero);
            obj.depths = obj.depths(nonzero);
            obj.rank = obj.rank(nonzero);
            %update E:
            nonzero = obj.S~=0;
            obj.E = zeros(obj.NR);
            obj.E(nonzero) = obj.S(nonzero)./obj.L(nonzero);
            
            obj.initVectorNormalization();
        end
        
        function idx = getSegmentsPixelIdxList(obj, s)
            idx = [];
            for i = 1:numel(s)
                ind = s(i);
                idx = [idx;...
                 cat(1,obj.s0.st(obj.hass(ind,1:obj.hasn(ind))).PixelIdxList)];
            end
        end
        
        %% Appearance-based similarity functions

        %ii is 1xN samples; returns NxD
        %ii should contain [i ineighbors]
        function V = getNormalizedVectors(obj, i0, ii)
            ni = size(ii,2);
            data = zeros(ni, 18); %7+7 + 4
            ii = ii';
            jj = i0*ones(ni,1);
            
            We = obj.AR([ii,jj]); %two cols for i,j
            sumWe = sum(We,2);
            We = We ./ [sumWe, sumWe]; % norm to sum 1 on rows
            % EX , ie, means. Not using depths here, only area weights.
            data(:,1:7)  = ( (We(:,1)*ones(1,7)) .* obj.EX(ii,:) ) ...
                         + ( (We(:,2)*ones(1,7)) .* obj.EX(jj,:) );
            newEX2       = (((We(:,1)*ones(1,7)) .* sqrt(obj.EX2(ii,:)))...
                         +  ((We(:,2)*ones(1,7)) .* sqrt(obj.EX2(jj,:)))...
                            ).^2;
            data(:,8:14) = real(sqrt(newEX2 - data(:,1:7).^2)); %STD
            % Centroids
            data(:,17:18)   = [... % weighted mean on x and on y
                  sum(We .* [obj.CC(ii,1), obj.CC(jj,1)], 2),... % CCx
                  sum(We .* [obj.CC(ii,2), obj.CC(jj,2)], 2)];   % CCy
            bbulx  = [obj.BB(ii,1) obj.BB(jj,1)];
            bbuly  = [obj.BB(ii,2) obj.BB(jj,2)];
            bbminx = min(bbulx,[],2); %col
            bbminy = min(bbuly,[],2); %col
            bbmaxx = max([bbulx(:,1)+obj.BB(ii,3), bbulx(:,2)+obj.BB(jj,3)],[],2);
            bbmaxy = max([bbuly(:,1)+obj.BB(ii,4), bbuly(:,2)+obj.BB(jj,4)],[],2);
            % BB = [bbminx, bbminy, bbmaxx-bbminx, bbmaxy-bbminy];
            data(:,15:16) = [bbmaxx-bbminx, bbmaxy-bbminy];

            % Normalization
            V = (data - ones(ni,1)*obj.means) ./ (ones(ni,1)*obj.stds);
        end
        function initVectorNormalization(obj)
            D = [obj.EX, sqrt(obj.EX2 - obj.EX.^2), obj.BB(:,3:4), obj.CC];
            obj.means = mean(D,1);
            obj.stds  = std( D,1);
            obj.stds(obj.stds==0) = 100; 
            % Anyway mean would be 0, then 0/100=0. Very infrequent
        end
        
        function precalcW(obj, k, mu, sigma)
            if exist('mu','var')
                obj.means = mu;
                obj.stds = sigma;
            end
            obj.W = -1*ones(obj.NR);
            
            %first get all 2nd neighbors if necessary
            A2 = obj.A;
            if k==2
                for i=1:obj.NR
                    nn = find(obj.A(i,:));
                    idx = any(obj.A(nn,:),1);
                    idx(nn) = true; idx(i) = false;
                    %nn = find(idx);d
                    A2(i,:) = idx;
                end
            end
            
            for i=1:obj.NR
                ineighbors = find(A2(i,:)); %ineighbors = obj.getKstepsNeighbours(i,k);
                V = obj.getNormalizedVectors(i, [i ineighbors]);
                V(:,8:14) = V(:,8:14)*0.75; %STD
                V(:,end-3:end-2) = V(:,end-3:end-2)*0.2; %W,H
                V(:,end-1:end  ) = V(:,end-1:end  )*0.2; %X,Y
                nn = size(V,1)-1;
                dif   = V(2:end,:) - ones(nn,1)*V(1,:);
                obj.W(i, ineighbors) ...
                      = sqrt(sum(dif.^2,2))' + obj.edgepriorW * ...
                        (~obj.A(i,ineighbors) + obj.E(i,ineighbors));
                % W is including edgeness.
                % to exclude edgeness: W2=obj.W-obj.edgepriorW*(~obj.A + obj.E)
            end

        end
                
        function [variances entropies nsuperpixels] ...
                 = getIntraSegmentVariances(obj, segments, meansvars)
            if ~exist('meansvars','var')
                meansvars = [obj.s0.EX, obj.s0.EX2 - obj.s0.EX .^2];
            end
            n = numel(segments);
            variances = zeros(n,14);
            entropies = zeros(1,n);
            nsuperpixels = zeros(1,n);
            for i=1:n
                s = segments(i);
                ns = obj.hasn(s);
                nsuperpixels(i) = ns;
                supps = obj.hass(s,1:ns);
                %variances(i,:) = var(meansvars(supps),[],2);
                if ns<2
                    entropies(i) = 0;
                    continue;
                end
                if ns>=5
                    ka = 5;
                else
                    ka = ns;
                end
                D = 7; 
                NNdist = dist(meansvars(supps,:)');
                H = -psi(ka) + psi(ns) + D/2*log(pi) - gammaln(D/2+1) + D/ns*sum(log(NNdist(NNdist~=0)));
                entropies(i) = H;
            end
        end
        
        
        function [sumw weights] = getFirstLevelWeightsOld(obj, s1, segment)
            ss0 = sort(obj.hass(segment,1:obj.hasn(segment)));
            ss0l = false(1,obj.s0.NR); ss0l(ss0) = true;
            weights = ones(obj.s0.NR)*NaN;
            
            for si = 1:numel(ss0) %for each s1 element in the segment,
                s = ss0(si);
                
                snn = ss0(logical(obj.s0.A(s,ss0l))); %the neighbors of the element
                for n=snn %and then the neighbors of the neighbors
                    snn = [snn ss0(logical(obj.s0.A(n,ss0l)))];
                end
                weights(s, snn) = s1.W(s,snn); % * correc;
            end
            weights =  (weights+weights')/2;
            
            sumw = sum(weights(~isnan(weights)))/sum(~isnan(weights(:)));
        end
        
        %In weights there are negative values which must be ignored
        function [sumw weights] = getFirstLevelWeights(obj, segment)
            ss0 = sort(obj.hass(segment,1:obj.hasn(segment)));
            ss0l = false(1,obj.s0.NR); ss0l(ss0) = true;
            %I should do this difference only once, outside this function:
            %It means I keep edgeness between 1st negihbors
            weights = obj.s0.s1.W(ss0l,ss0l) - obj.edgepriorW*(~obj.s0.s1.A(ss0l,ss0l));% + obj.s0.s1.E(ss0l,ss0l));
            weights =  (weights+weights')/2;
            weights(logical(eye(size(weights)))) = 0;
            sumw = sum(weights(weights>=0))/sum(weights(:)>=0);
        end
        
        %Heterogeneity between first and second NN only
        %In weights there are negative values which must be ignored
        function [htg] = heterogeneity(obj, segments, percentile)
            if isempty(obj.heterogen)
                obj.heterogen = -1*ones(1,obj.NR); %temporaries
            end
            for segment = segments(obj.heterogen(segments)==-1)
                
                [s, wei] = obj.getFirstLevelWeights(segment);
                obj.heterogen(segment) = prctile(wei(wei>=0),percentile);
            end            
            htg = obj.heterogen(segments);
        end
        
        %Heterogeneity between all superpixels
        % segments are indexes of segments at the current level
        % meansvars are m x 16, for all the m superpixels at level 0,
        % containing at each row the 7 means and 7 variances from EX,EX2
        function [htg] = heterogeneity2(obj, segments, percentile, meansvars)
            n = numel(segments);
            htg = zeros(1,n);
            for i=1:n
                s = segments(i);
                ns = obj.hasn(s);
                supps = obj.hass(s,1:ns);
                wei = dist(meansvars(supps,:)');
                htg(i) = prctile(wei(wei>=0),percentile); 
            end
        end
        
        function [entropies] = get2NNentropies(obj, segments, meansvars)
            if ~exist('meansvars','var')
                meansvars = [obj.s0.EX, obj.s0.EX2 - obj.s0.EX .^2, obj.s0.CC];
            end
            extreme = 1000000;
            n = numel(segments);
            entropies = zeros(1,n);
            for i=1:n
                s = segments(i);
                ns = obj.hasn(s);
                ss0 = sort(obj.hass(s,1:ns));
                ss0l = false(1,obj.s0.NR); ss0l(ss0) = true;
                weights = obj.s0.s1.W(ss0l,ss0l) - obj.edgepriorW*(~obj.s0.s1.A(ss0l,ss0l));% + obj.s0.s1.E(ss0l,ss0l));
                maxka = sum(weights>0, 2);
                weights(logical(eye(size(weights)))) = 0;
                weights(weights<0) = extreme;
                NNdistall = sort(weights,2); %Idea: sort by 1st or 2nd NN first, then distance.
                D=16;
                if ns<2
                    continue;
                end
                mka = max(maxka);
                for ka = 1:mka
                    NNdist = NNdistall(:,ka+1);
                    NNdist = NNdist(NNdist~=extreme);
                    ns = numel(NNdist);
                    H = -psi(ka) + psi(ns) + D/2*log(pi) - gammaln(D/2+1) + D/ns*sum(log(NNdist(NNdist~=0)));
                    entropies(i) = entropies(i) + H;
                end
                entropies(i) = entropies(i) / mka;
                
            end
        end
        
        function COV = getCovarAppCoord(obj, s)
            ih = obj.s0.ih;
            nr = obj.hasn(s);
            idxs = [];
            for i = 1:nr
                idxs = [idxs; obj.s0.st(obj.hass(s,i)).PixelIdxList(:) ];
            end
            [I,J] = ind2sub(obj.s0.RC,idxs);
            % Nsamples x Ndimensions
            D = [ih.labl(idxs),...
                 ih.laba(idxs),...
                 ih.labb(idxs),...
                 ih.gx(idxs)*0.002,...
                 ih.gy(idxs)*0.002,...
                 I*0.006,...
                 J*0.006
                 %ih.gxx(idxs),...
                 %ih.gyy(idxs),...
                 ];
            COV = reshape(cov(D),1,[]);
        end
        
        %% Plot functions
        function segments = getSegmentsPlot(obj)
            segments = zeros(obj.s0.RC);
            
            for s=1:obj.NR
                segments(cat(1,obj.s0.st(obj.hass(s,1:obj.hasn(s))).PixelIdxList)) = s;%rand(1);
            end
        end
        
        function segments = getSegmentsInteriorMask(obj, thickness)
            ofs = thickness; % 2 or 4
            a = obj.getSegmentsPlot();
            v = a(:,1:end-ofs) == a(:,ofs+1:end);
            h = a(1:end-ofs,:) == a(ofs+1:end,:);
            segments = [false(ofs/2,obj.s0.RC(2)-ofs); (v(ofs:end-1,:) & h(:,ofs:end-1)); false(ofs/2,obj.s0.RC(2)-ofs)];
            segments = [false(ofs/2,obj.s0.RC(1)); segments'; false(ofs/2,obj.s0.RC(1))]';
        end
           
        
        function disp = plotGraphsHeat(obj,H)
            mui = H;
            s = sort(H);
            mindif = min(abs(s(2:end)-s(1:end-1)));
            noise = rand(size(H))*(mindif*0.2);
            mui = mui + noise;
            disp = zeros(obj.s0.RC);
            for ind=1:obj.NR
                if mui(ind) ~= 0
                disp(cat(1,obj.s0.st(obj.hass(ind,1:obj.hasn(ind))).PixelIdxList)) = mui(ind);
                end
            end
        end
        
        function disp = plotGraphsHeatOriginal(obj,H)
            mui = H;
            s = sort(H);
            disp = zeros(obj.s0.RC);
            for ind=1:obj.NR
                if mui(ind) ~= 0
                disp(cat(1,obj.s0.st(obj.hass(ind,1:obj.hasn(ind))).PixelIdxList)) = mui(ind);
                end
            end
        end
        
        
        % level=45; region=1; imshow(segh(level).markRegions(setdiff(1:segh(level).NR, region), imread(imagepath), 255,255,255));
        function out = markRegions(obj, regionidxlist, img, R, G, B)
            %if ~exist('B','var'); R = 0; G = 0; B = 255; end
            if ~exist('img','var'); img = zeros([obj.s0.RC 3]); end
            out = img;
            [rows, cols, ~] = size(img);
            ind = obj.getSegmentsPixelIdxList(regionidxlist);
            
            %inoise = rand(size(ind))<0.2;
            %vnoise = uint8(rand(1,sum(inoise))*255);
            %indnoise = ind(inoise);
            
            out(ind) = out(ind)*5;
            out(ind+numel(img(:,:,1))) = out(ind+numel(img(:,:,1)))*0.7;
            out(ind+2*numel(img(:,:,1))) = out(ind+2*numel(img(:,:,1)))*0.7;
            
            out(ind) = R;
            out(ind+numel(img(:,:,1))) = G;
            out(ind+2*numel(img(:,:,1))) = B;
            
        end
        
        function out = maskRegions(obj, regionidxlist)
            out = logical(zeros(obj.s0.RC));
            ind = obj.getSegmentsPixelIdxList(regionidxlist);
            out(ind) = 1;
        end
        
    end
end
