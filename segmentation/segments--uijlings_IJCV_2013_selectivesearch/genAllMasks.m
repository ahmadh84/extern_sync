%Takes selective search's output from selective_search_segs and generates
% object segmentations
%
% @authors:     Ahmad Humayun
% @contact:     ahumayun@cc.gatech.edu
% @affiliation: Georgia Institute of Technology
% @date:        Fall 2014
function [ masks ] = genAllMasks(segs)
    num_cfgs = numel(segs.hierarchies);
    
    testBlobs = cell(num_cfgs,1);
    % iterate over all the configurations, and generate blobs for them
    for idx = 1:num_cfgs
        testBlobsT = cell(length(segs.hierarchies{idx})+1, 1);
        
        % get the initial blobs (superpixels) and the first hierarcy blobs
        [~, testBlobsT{end}, testBlobsT{1}] = ...
            RecreateBlobHierarchyIndIm(segs.sp_maps(:,:,idx), ...
                                       segs.sp_blob_boxes{idx}, ...
                                       segs.hierarchies{idx}{1});
        % get the blobs from all the other hierarchies
        for j = 2:length(segs.hierarchies{idx}) % Without initial blobs here
            [~, ~, testBlobsT{j}] = ...
                RecreateBlobHierarchyIndIm(segs.sp_maps(:,:,idx), ...
                                           segs.sp_blob_boxes{idx}, ...
                                           segs.hierarchies{idx}{j});
        end
        % merge all the blobs for this configuration
        testBlobs{idx} = cat(1, testBlobsT{:});
%         fprintf('%d\n', numel(testBlobs{idx}));
    end
    % merge all blobs across all configurations
    testBlobs = cat(1, testBlobs{:});
    
    % reorder and filter blobs according to selective_search_segs
    testBlobs = testBlobs(segs.order_idxs);
    
    % sanity check
    b = vertcat(testBlobs{:}); b = vertcat(b.rect);
    assert(all(all(segs.bb_boxes == b)));
    
    sz = size(segs.sp_maps);
    masks = false(sz(1), sz(2), numel(testBlobs));
    % now generate mast for each blob
    for m_idx = 1:numel(testBlobs)
        blob = testBlobs{m_idx};
        masks(blob.rect(1):blob.rect(3), ...
              blob.rect(2):blob.rect(4), m_idx) = blob.mask;
    end
end

