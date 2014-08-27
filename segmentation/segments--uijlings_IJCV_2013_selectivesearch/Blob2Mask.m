function [ mask, imout ] = Blob2Mask(blob, im)
    mask = false(size(im,1), size(im,2));
    mask(blob.rect(1):blob.rect(3), blob.rect(2):blob.rect(4)) = blob.mask;
    
    if nargout > 1
        imout = im;
        % set pixels not in mask to zero
        imout(repmat(~mask, [1 1 size(im,3)])) = 0;
        imout = imout(blob.rect(1):blob.rect(3), blob.rect(2):blob.rect(4), :);
    end
end

