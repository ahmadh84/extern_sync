function [ masks ] = savecands2masks( candidates )
%SAVECANDS2MASKS Summary of this function goes here
%   Detailed explanation goes here

    n_cands = numel(candidates.labels);
    sz = size(candidates.superpixels);
    masks = false(sz(1),sz(2),n_cands );
    for i = 1:n_cands
        masks(:,:,i) = ismember(candidates.superpixels, candidates.labels{i});
    end
end

