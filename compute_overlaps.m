function f = compute_overlaps(masks, GT, care_part, type)
% Can compute different error measures given an output mask and a GT. It
% can also compute errors for multiple masks at the same time. The output
% is a vector of error scores of length N where N is the number of masks 
% provided.
%
% @args:
%   masks: is a logical array (size of the image) which gives the mask
%      output by an algorithm, whose overlap needs to be computed with the
%      provided ground-truth. You can provide multiple masks by
%      concatenating them on the 3rd dimension
%
%   GT: is an array (size of the image) giving the ground-truth over which 
%       all the masks are scored - which are 0 for negatives and 1 for 
%       positives
%
%   care_part: is a logical array (size of the image) telling the parts 
%       which you want to count toward TN/FP/FN/TP.
%
%   type (<'overlap'>, 'intersect_gt', 'centroid_displacement'): specifies
%       the type of overlap you want to compute. 'overlap' computes
%       TP/(FP+FN+TP). 'intersect_gt' computes the precision TP/(FP+TP)

    gtcen = regionprops(GT,'Centroid');
    f = zeros(size(masks,3),1);
    % This needs memory though, decide later if it's worth it
%    masks = uint8(masks);

    % this operation would set all masks to this: 0=TN, 1=FP, FN=2, TP=3
    % (all others values are don't care)
    masks_pl_gt = bsxfun(@plus, uint8(GT) * 2, uint8(masks));
    masks_pl_gt = reshape(masks_pl_gt, size(masks_pl_gt,1) * size(masks_pl_gt,2), size(masks_pl_gt,3));
    
    if exist('care_part','var') && ~isempty(care_part)
        all_is = histc(masks_pl_gt(care_part,:), 0:3, 1);
    else
        all_is = histc(masks_pl_gt, 0:3, 1);
    end
    all_is = all_is';
    if ~exist('type','var') || isempty(type) || strcmp(type,'overlap')
        f = all_is(:,4) ./ (sum(all_is(:,2:4),2) + eps);
    elseif strcmp(type,'intersect_gt')
        f = all_is(:,4) ./ (all_is(:,2) + all_is(:,4) + eps);
    elseif strcmp(type,'centroid_displacement')
        for i=1:size(masks,3)
            maskcen = regionprops(masks(:,:,i),'Centroid');
            gtcen = gtcen.Centroid;
            maskcen = maskcen.Centroid;
            f(i) = norm(gtcen - maskcen);
        end
    end
end