% Canidate Regions - version 1.0.  Code distributed under LGPL. 
% August 2014, University of California, Los Angeles
%
% Please, cite this work if you use the code: 
%   B. Bonev, A.L. Yuille, "A Fast and Simple Algorithm for Producing Candidate Regions", 
%   ECCV 2014, Zürich, Switzerland, 6th-12th September, 2014
%

function demo()
basefolder = [pwd() '/../']; %if you use a relative path it doesn work for opts.modelDir
%% External libraries
%You may need to modify the folders here, if you use your own libraries
%instead of the ones that can be downloaded and placed in the lib/ folder.
%VL feat library
addpath('/home/ahumayun/extern_src/toolboxes/vlfeat/toolbox/');
vl_setup();
%Piotr Toolbox library
addpath(genpath('/home/ahumayun/extern_src/toolboxes/piotr_toolbox/'));
%Dollar's ICCV 2013 Structured Edge Detector code initialization:
addpath([basefolder '/lib/piotr/release/']);
currentfolder = cd();
%Load Dollar's sketch token model:
opts=edgesTrain();                % default options (good settings)
opts.modelDir=[basefolder '/lib/piotr/release/models/'];% model in models/forest
opts.modelFnm='modelFinal';       % model name
opts.nPos=5e5; opts.nNeg=5e5;     % decrease to speedup training
opts.useParfor=0;                 % parallelize if sufficient memory
edgesmodel=edgesTrain(opts);
edgesmodel.opts.multiscale=0;     % for top accuracy set multiscale=1
edgesmodel.opts.nTreesEval=4;     % for top speed set nTreesEval=1
edgesmodel.opts.nThreads=4;       % max number threads for evaluation
edgesmodel.opts.nms=0;            % set to true to enable nms (fairly slow)
cd(currentfolder);

%% Input and output folders
resultsfolder = [basefolder '/output/'] 
resultsfoldercr = [resultsfolder '/candidate_regions/'] 
resultsfolderbb = [resultsfolder '/boxes/'] 
imfolder = '/data/images/everingham_IJCV_2010_pascalvoc/JPEGImages/';
imageset = '/data/algo_results/segmentation/segments--bonev_ECCV_2014_candidateregions/input/ImageSet/sample.txt';
%addpath([basefolder 'Datasets/Pascal/VOCdevkit/VOCcode']); %For groundtruth
%VOCinit2011;
groundtruthfolder = '/data/algo_results/segmentation/segments--bonev_ECCV_2014_candidateregions/input/Groundtruth/class/';
groundtruthobjectfolder = '/data/algo_results/segmentation/segments--bonev_ECCV_2014_candidateregions/input/Groundtruth/instance/';

%% Some parameters
mmasks = true; % Produces Candidate-Regions as described in the paper
bboxes = true; % Produces bounding boxes, optional
evaluateWithGroundtruth = true; % Finds the best segments and their Intersection over Union (IoU), included in the filenames 

% Segmentation: 
slicRegu = 0.1;   % smaller -> more jagged SLIC segments
ratioMerge = 0.3; % Amount of segments which will be selected for merging
                  % at each level, from 0 to 1. Larger -> less levels
nfewsegments = 50;% At how many segments to stop using k=2 and start k=1 (for merges)
edgesweight =1.35;% How important the edgemap will be in the dissimilarity measure


%% Other initializations
format short g;
clear evaluateBoxes; %these functions have persistent variables
clear evaluateMasks;
evaluateBoxes('clear');
evaluateMasks('clear');
trn_file = fopen(imageset);
files1 = textscan(trn_file,'%s','Delimiter','\n');
files = files1{1};

average_totaltime = 0;
average_number_candregions_andboxes = [0 0];
ntimes = 0; 

classnames = {'plane','bicycle','bird','boat','bottle','bus','car','cat',...
'chair','cow','d. table','dog','horse','motorbike','person',...
'pot. plant','sheep','sofa','train','tv',...
'bag','bed','bench','book','building','cabinet','ceiling','clothes',...
'pc','cup','door','fence','floor','flower','food','grass',...
'ground','keyboard','light','mountain','mouse','curtain','platform',...
'sign','plate','road','rock','shelves','sidewalk','sky','snow',...
'table','track','tree','truck','wall','water','window','wood'};

for i = 1:numel(files)
    
%% Load groundtruth
    
    sprintf('%d of %d:',i,numel(files))
    
    name = files{i}
    imgn = name(1:11);
    imagepath = [imfolder imgn '.jpg'];

    if evaluateWithGroundtruth
        try
            doEvaluate = true;
            gtruthclass = imread([groundtruthfolder imgn '.png']);
            try
                gtruthobject = imread([groundtruthobjectfolder imgn '.png']);
            catch excep
                gtruthobject = zeros(size(gtruthclass));
                display('Warning: missing "instances" groundtruth for this image')
            end
            segt = superpixelsHandle();
            [segments, correspondences] = segt.renumberSegments(gtruthclass+1, gtruthobject);
            correspondences = correspondences - 1;
            segt.setSegments(segments);
            segt.regionprops();
            segt.adjacency(zeros(segt.RC));
        catch err
            doEvaluate = false;
            disp(err);
            disp(err.message);        
        end
    else
        doEvaluate = false;
    end
%% Segmentation hierarchy
    originalimage = imread(imagepath);
    tic
    slicSize = floor(size(originalimage,1)*size(originalimage,2)*0.00011); if slicSize<10, slicSize=10;end
    [segh] = hdsegments3(originalimage, false, edgesmodel, edgesweight, slicSize, slicRegu, ratioMerge, nfewsegments);
    segmentation_time = toc
        
%% Segments selection by entropy criterion
    tic
    minHthres = 25;
    [nselected segh] = selectionByEntropy(segh,minHthres);
    selection_time = toc    
            
%% Combinations of up to 3 selected-segments

    tic
    nlevels = numel(segh);
    % supersegments generation:    
    combs3 = supersegmentCombinations(segh, nlevels-10:nlevels-1, 5, 3, 10000, 0.15); 
    combs2 = supersegmentCombinations(segh, 10:nlevels-1, 5, 2); 
    combs1 = supersegmentCombinations(segh, 1:9, 0, 1);      
    combs = {combs1{:},combs2{:},combs3{:}};        

    %beforefilter = numel(combs)
    [combs, combx] = filterCombs(combs, segh, 0.05, mmasks, bboxes);
    nsegs_nbox = [size(combs,2), size(combx,2)];

    %reduce if too many parts
    if mmasks && nsegs_nbox(1) > 1200
        combs3 = supersegmentCombinations(segh, nlevels-15:nlevels-1, 15, 3, 30000, 0.2); 
        combs2 = supersegmentCombinations(segh, 21:nlevels-1, 10, 2); 
        combs1 = supersegmentCombinations(segh, 1:20, 0, 1);      
        combs = {combs1{:},combs2{:},combs3{:}};
        [combs, combx] = filterCombs(combs, segh, 0.05, mmasks, bboxes);
        nsegs_nbox = [size(combs,2), size(combx,2)];
    end
     
%% Save results and evaluate
    combinations_time = toc


    % Save and evaluate Boxes
    if bboxes
        [parts, partsip, partslevel] = saveResultingBoxes(segh,combx, resultsfolderbb,imgn);
        if doEvaluate
            evaluateBoxes(i, parts, partsip, partslevel, segt, correspondences, originalimage, [resultsfolderbb '/eval/'], imgn, '', classnames)
        end
    end
    % Save and evaluate Masks
    if mmasks
        [parts, partsip, partslevel] = saveResultingMasks(segh,combs, resultsfoldercr,imgn);
        if doEvaluate
            evaluateMasks(i, parts, partsip, partslevel, segt, correspondences, originalimage, [resultsfoldercr '/eval/'], imgn, '', classnames)
        end
    end

    elapsedtime = segmentation_time + selection_time + selection_time;
    average_totaltime = (elapsedtime + average_totaltime*ntimes)/(ntimes+1)
    average_number_candregions_andboxes = (nsegs_nbox + average_number_candregions_andboxes*ntimes)/(ntimes+1)
    ntimes = ntimes + 1;

end %for

display('Results are stored in the output folder: ');
display('The .mat files contain all the masks for the CR (candidate regions).');
display('The eval subfolder shows only the best segments after evaluation (if available).'); 
display('To check if the produced results are the expected, compare with the output_expected folder provided with this demo.');
display('Average totaltime should be ~4.0s, using Linux 64bits, Matlab2013b on a laptop with Intel(R) Core(TM) i7-2640M CPU @ 2.80GHz.');
display('Average #candidate regions and #boxes should be 481 and 311 for the 7 sample images provided with this demo.');

end %run
%% Segmentation evaluation functions

function [IoU pixI pixU] = evaluateSupersegment(segmask, segt)
    %This is IoU of MASKS evaluation.
    % For each s in segt use its associated segments in seg to get the
    % union and intersection
    IoU = zeros(1,segt.NR);
    pixI= zeros(1,segt.NR);
    pixU= zeros(1,segt.NR);
    for isegt = 1:segt.NR
        segtmask = false(size(segmask));
        segtmask( segt.st(isegt).PixelIdxList(:)' ) = true;
        pixI(isegt) = sum(sum(segtmask & segmask));
        pixU(isegt) = sum(sum(segtmask | segmask));
        IoU(isegt) =  pixI(isegt) / pixU(isegt);
        if pixI(isegt)==0
            pixU(isegt) = sum(sum(segtmask));
        end
    end
end
function [IoU pixI pixU] = evaluateSupersegmentBBox2(bbxs, segt)
    %This is IoU of MASKS evaluation.
    % For each s in segt use its associated segments in seg to get the
    % union and intersection
    IoU = zeros(1,segt.NR);
    pixI= zeros(1,segt.NR);
    pixU= zeros(1,segt.NR);
    for isegt = 1:segt.NR
        segtmask = segt.segments == isegt;
        anyy = any(segtmask,2); anyx = any(segtmask,1);
        bbxt = [find(anyy,1,'first'), find(anyy,1,'last'), ...
               find(anyx,1,'first'), find(anyx,1,'last')];
        %bbxt = [miny maxy minx maxx];   
        
        bbxt = [bbxt(1),bbxt(3), bbxt(2)-bbxt(1)+1, bbxt(4)-bbxt(3)+1];
        
        pixI(isegt) = rectint(bbxt,bbxs);
        pixU(isegt) = bbxt(3) * bbxt(4) + bbxs(3) * bbxs(4) - pixI(isegt);
        
        IoU(isegt) =  pixI(isegt) / pixU(isegt);
        if pixI(isegt)==0
            pixU(isegt) = sum(sum(segtmask));
        end
    end
end

%% Supersegment combinations
function listall = supersegmentCombinations(segh, levels, maxdepth, maxlen, minsize, minrelative)
    if ~exist('minsize', 'var')
        minsize = 0;
    end
    if ~exist('minrelative', 'var')
        minrelative = 0;
    end
    listall = {};
    for level = levels
      for s = segh(level).parts
        ssize = segh(level).AR(s);
        if maxlen>1 & ssize >= minsize
            [nns nnlevels] = segh(level).getMultilevelNN(s, maxdepth, 0.3);
            %filter out combs with too small neighbours
            keep = false(size(nns));
            for i=1:numel(nns)
                if segh(nnlevels(i)).AR(nns(i)) >= minrelative*ssize
                    keep(i) = true;
                end
            end
            nns = nns(keep);
            nnlevels = nnlevels(keep);
        else
            nns = [];
        end
        list = cell(1); 
        list{1} = [s level];
        if ~isempty(nns)
            %Get a list of combinations for segment s with its neighbours
            for i=1:size(nns,2)
                for l = 1:size(list,2)
                    listl = list{l};
                    if size(listl,1)<maxlen
                      list{end+1} = [listl; nns(i) nnlevels(i)];
                    end
                end
            end
        end
        for l = 1:size(list,2)
            si = size(list{l},2);
            if ~isempty(list{l})
                listall = {listall{:}, list{l}};
            end
        end
      end        
    end
end

function [outcombsmasks, outcombsboxes] = filterCombs(combs, segh, percentage, withMasks, withBboxes)
    outcombsmasks = {};
    outcombsboxes = {};
    nl = numel(segh);
    nc = numel(combs);
    masks = false(nc, segh(1).s0.NR); %Superpixel representation
    for i=1:nc
        c = combs{i};
        for j=1:size(c,1)
            l = c(j,2);
            s = c(j,1);
            masks(i,segh(l).hass(s, 1:segh(l).hasn(s) )) = 1;
        end
    end
    if withBboxes
        bbxs = zeros(nc,4);
        for i=1:nc
            bxslic = cat(1, segh(1).s0.st(masks(i,:)).BoundingBox);
            yx1 = min(bxslic(:,1:2),[],1);
            yx2 = max(bxslic(:,1:2),[],1); %should add the w and h of this one; for checking size doesn't matter
            bbxs(i,:) = [yx1, yx2-yx1];
        end
        [uni idx] = unique(bbxs,'rows');
        outcombsboxes = {combs{idx}};
    end
    if withMasks
        if percentage == 0
            [uni idx] = unique(masks,'rows');
            outcombsmasks = {combs{idx}};
        else %That's a lot slower:
            idxRepeatedMasks = false(nc,1);
            for i=nc:-1:1
                if idxRepeatedMasks(i)
                    continue
                end
                p = (masks(i,:));
                sizep = sum(p(:));
                xors = xor(masks(1:i-1,:), repmat(p, [i-1 1]));
                sizes = sum(xors,2); % 0 means identical

                idxRepeatedMasks(1:i-1) = idxRepeatedMasks(1:i-1) | sizes < (percentage)*sizep;
            end
            outcombsmasks = {combs{~idxRepeatedMasks}};
        end
    end
end

%% Segments selection

function [totalselected segh] = selectionByEntropy(segh, minHgain)
        %Here the segments to evaluate are those in segh(:).parts
        oldentropy = zeros(segh(1).NR,1);
        oldnsupp = ones(segh(1).NR,1);
        totalselected = 0;
        % These mean and std values are calculated from a random set of 1000 images
        % Contents of the vectors mm and std: 
        % 7 means (l, a, b, gx, gy, gxx, gyy), 7 vars (same), 1 w, 1 h, 1 cx, 1 cy  
        mm = [0.351532485888255,0.333408080050130,0.490875888844171,-0.0491434962780057,-0.111764742975605,-9.24717475374528,-4.52616528950800,0.0237856052152852,0.0102335143848912,0.00444232288707433,10983.0794976152,6572.79632471097,281046.246628482,180480.202872316, 247.645146337971,181.435789545103];
        std = [0.158234017317445,0.170716613550236,0.203984850279852,19.4411782682435,12.1987204912077,311.238338932908,173.139905539515,0.0329851872653395,0.0170696113141047,0.00580531768626500,18362.7598836125,12424.3232121515,417480.731856197,328792.460145005, 141.824505557249,103.520989143669];
        segh(1).precalcW(2,[mm, 255, 170], [std, 146, 96]);
        meansvars = [segh(1).s0.EX, segh(1).s0.EX2 - segh(1).s0.EX.^2 ];
        meansvars = (meansvars-repmat(mm(1:14),size(meansvars,1),1))./repmat(std(1:14),size(meansvars,1),1);
        for level = 2:numel(segh)
            %level
            seg = segh(level);
            
            [~, entropy, nsupp] = seg.getIntraSegmentVariances(1:seg.NR,meansvars);

            merged = find(seg.hads(:,2)>0); 
            large = oldnsupp(seg.hads(merged,1))>1 & oldnsupp(seg.hads(merged,2))>1;
            ml = merged(large);
            badmerged = ~diag(segh(level-1).A(seg.hads(merged,1),seg.hads(merged,2)))';
            if ~isempty(merged)
                if ~isempty(ml)
                    meanentropies = sum(oldentropy(seg.hads(ml, :))   ,2);
                    difs = (entropy(ml) - meanentropies'); 
                    badmerged(large) = badmerged(large) | difs>minHgain;
                end

                %goodsegments and seg.hads are in level-1
                %merged(badmerged) is in level-0
                goodsegments = reshape(seg.hads(merged(badmerged),:),1,[]);
                totalselected = totalselected + numel(goodsegments);

                segh(level-1).parts = goodsegments;
            else
                segh(level-1).parts = [];
            end
            
            oldentropy = entropy;
            oldnsupp = nsupp;
        end
end

%% Results saving functions

function [parts, partsip, partslevel] = saveResultingMasks(segh,combs, resultsfolder, imgn, writeEachSegment)
    if ~exist('writeEachSegment','var')
        writeEachSegment = false;
    end
    boxsuffix = '';
    totalparts = numel(combs);
    partsip = zeros(2,totalparts);
    partslevel = zeros(2,totalparts);
    parts = false([totalparts, segh(1).s0.RC]); %logical zeros
    % Data structures for the combinations
    for ic=1:size(combs,2)
        if isempty(combs{ic})
            continue
        end
        combsic = combs{ic};
        combsn = combsic(:,1);
        combsl = combsic(:,2);

        part = false(segh(1).s0.RC);
        for iii = 1: size(combsn,1)
            part(segh(combsl(iii)).getSegmentsPixelIdxList(combsn(iii))')...
                = true;
        end
        %disksize=round(sqrt(sum(sum(part)))/40);
        %disk = strel('disk',disksize);
        %parts(ic,:,:) = imerode(fillgaps(imdilate(part,disk), 2000, 0.35), disk);
        parts(ic,:,:) = part; %fillgaps(part, 2000, 0.35);

        [~,a] = sort(combsl,'descend');
        na = numel(a);
        partslevel(1:na,ic) = combsl(a);
        partsip(1:na,ic) = combsn(a);
        
        if writeEachSegment %(!)
         imwrite(part,[resultsfolder imgn sprintf('%s_L%02d_P%02d_L%02d_P%02d',...
            '', partslevel(1,ic), partsip(1,ic), partslevel(2,ic), partsip(2,ic)) '.png'],'png');
        end
    end
    save([resultsfolder imgn boxsuffix '.mat'],'parts');

end
function [parts, partsip, partslevel] = saveResultingBoxes(segh,combs, resultsfolder, imgn)
    boxsuffix = 'BB';
    totalparts = numel(combs);
    parts = zeros(totalparts, 4);
    % Data structures for the combinations
    for ic=1:size(combs,2)
        if isempty(combs{ic})
            continue
        end
        combsic = combs{ic};
        combsn = combsic(:,1);
        combsl = combsic(:,2);

        part = false(segh(1).s0.RC);
        for iii = 1: size(combsn,1)
            part(segh(combsl(iii)).getSegmentsPixelIdxList(combsn(iii))')...
                = true;
        end
        anyy = any(part,2); anyx = any(part,1);
        bbx = [find(anyy,1,'first'), find(anyy,1,'last'), ...
               find(anyx,1,'first'), find(anyx,1,'last')];
        bbx = [bbx(1),bbx(3), bbx(2)-bbx(1)+1, bbx(4)-bbx(3)+1];
        parts(ic,:) =  bbx;

        [~,a] = sort(combsl,'descend');
        na = numel(a);
        partslevel(1:na,ic) = combsl(a);
        partsip(1:na,ic) = combsn(a);
    end
    save([resultsfolder imgn boxsuffix '.mat'],'parts');

end

function evaluateMasks(i, parts, partsip, partslevel, segt, correspondences, originalimage, resultsfolder, imgn, experimentsuffix, classnames)
    persistent tableS;
    persistent fillratio;
    persistent minfillratio;
    persistent ntimes;
    if nargin==1
        fillratio = 0;
        minfillratio = 1;
        ntimes = 0;
        tableS = [];
        fprintf('Clearing tableS for evaluateMasks\n');
        return
    end
    
    if size(partslevel,1)==1
        partslevel = [partslevel; zeros(2,size(partslevel,2))];
        %partsip = [partsip; partsip];
    end
    if size(partslevel,1)==2
        partslevel = [partslevel; zeros(1,size(partslevel,2))];
        %partsip = [partsip; partsip];
    end
    if isempty(tableS)
        tableS = [];
    end
    % Evaluate combinations
    IoU = zeros(size(parts,1), size(correspondences,2));
    pixi= zeros(size(parts,1), size(correspondences,2));
    pixu= zeros(size(parts,1), size(correspondences,2));
    for ic=1:size(parts,1)
        [ious pixis pixus] = evaluateSupersegment(squeeze(parts(ic,:,:)), segt);
        IoU(ic,:) = ious;
        pixi(ic,:) = pixis;
        pixu(ic,:) = pixus;
    end
    [iou, ibest] = max(IoU);
    if size(parts,1) == 0
        fprintf('Zero parts!');
        iou = [];
    end
    for j = 1:size(iou,2)
        object = correspondences(j);
        oarea  = segt.st(j).Area;
        bbarea = prod(segt.st(j).BoundingBox(3:4));
        %   1        2        3            4          5    
        %img_id, img_year, img_num, selected level, object,
        %
        %         6          7    8         9      10       11        12   
        % object area (gt), iou, nrs. combs_obj combs_img. level  bb area(gt) .
        %
        %      13    14      
        % pixelsI, pixelsU 
        tableS = [tableS; double(i), str2double(imgn(1:4)), str2double(imgn(6:end)), ...
                 double(-1), double(object), double(bbarea), ...
                 double(iou(j)), double(size(find(partslevel(:,ibest(j))~=0),1)), ...
                 double(size(find(IoU(:,j)),1)), double(size(parts,1)), ...
                 double(partslevel(1,ibest(j))), double(oarea), ...
                 double(pixi(ibest(j),j)), double(pixu(ibest(j),j))];
        if object>0 %&& object <=20
            bestL = partslevel(1,ibest(j)); bestip = partsip(1,ibest(j));
            %imwrite(im2bw(segh(bestL).markRegions(segh(bestL).parts(bestip)),0.1), ...

            part = squeeze(parts(ibest(j),:,:));
                        
            imwrite(part,...
             [resultsfolder imgn sprintf('_%s_id%02d_IoU%0.2f',...
             classnames{object},j,iou(j)), '.png'],'png');
        end
    end
    %save([resultsfolder 'tableS' experimentsuffix '.mat'],'tableS');
end
function evaluateBoxes(i, parts, partsip, partslevel, segt, correspondences, originalimage, resultsfolder, imgn, experimentsuffix, classnames)
    persistent tableS;
    persistent colors;

    if nargin==1
        tableS = [];
        fprintf('Clearing tableS for evaluateBoxes\n');
        return
    end
    
    if size(partslevel,1)==1
        partslevel = [partslevel; partslevel];
        partsip = [partsip; partsip];
    end
    
    if isempty(tableS)
        tableS = [];
        colors = uint8(hsv*255);
    end
    boxsuffix = 'BB';
    % Evaluate combinations
    IoU = zeros(size(parts,1), size(correspondences,2));
    pixi= zeros(size(parts,1), size(correspondences,2));
    pixu= zeros(size(parts,1), size(correspondences,2));
    imboxes = originalimage;
    for ic=1:size(parts,1)
        [ious pixis pixus] = evaluateSupersegmentBBox2(parts(ic,1:4), segt);
        IoU(ic,:) = ious;
        pixi(ic,:) = pixis;
        pixu(ic,:) = pixus;
        y1 = parts(ic,1);
        x1 = parts(ic,2);
        y2 = sum([parts(ic,[1,3]), -1]);
        x2 = sum([parts(ic,[2,4]), -1]);
        color = colors(ceil(rand(1)*size(colors,1)),:);
        imboxes(y1:y2,[x1,x2],[1]) = color(1); 
        imboxes(y1:y2,[x1,x2],[2]) = color(2); 
        imboxes(y1:y2,[x1,x2],[3]) = color(3); 
        imboxes([y1,y2],x1:x2,[1]) = color(1); 
        imboxes([y1,y2],x1:x2,[2]) = color(2); 
        imboxes([y1,y2],x1:x2,[3]) = color(3); 
    end
    [iou, ibest] = max(IoU);
    if size(parts,1) == 0
        fprintf('Zero parts!');
        iou = [];
    end
    for j = 1:size(iou,2)
        object = correspondences(j);
        oarea  = segt.st(j).Area;
        bbarea = prod(segt.st(j).BoundingBox(3:4));
        %   1        2        3            4          5    
        %img_id, img_year, img_num, selected level, object,
        %
        %         6          7    8         9      10       11        12   
        % object area (gt), iou, nrs. combs_obj combs_img. level  bb area(gt) .
        %
        %      13    14      
        % pixelsI, pixelsU 
        tableS = [tableS; double(i), str2double(imgn(1:4)), str2double(imgn(6:end)), ...
                 double(-1), double(object), double(bbarea), ...
                 double(iou(j)), double(size(find(partslevel(:,ibest(j))~=0),1)), ...
                 double(size(find(IoU(:,j)),1)), double(size(parts,1)), ...
                 double(partslevel(1,ibest(j))), double(oarea), ...
                 double(pixi(ibest(j),j)), double(pixu(ibest(j),j))];
        if object>0 %&& object <=20
            part = parts(ibest(j),:);
            bbx  = [part(1:2), part(1:2)+part(3:4)-1];
            part = uint8(zeros(size(originalimage)));
            part(bbx(1):bbx(3),bbx(2):bbx(4),:) = uint8(originalimage(bbx(1):bbx(3),bbx(2):bbx(4),:));
            imwrite(part,...
             [resultsfolder imgn boxsuffix sprintf('_%s_id%02d_IoU%0.2f',...
             classnames{object}, j, iou(j)), '.png'],'png');
        end
    end
    %save([resultsfolder 'table' boxsuffix experimentsuffix '.mat'],'tableS');
    %imwrite(imboxes,...
    % [resultsfolder imgn boxsuffix '.png'],'png');
end









