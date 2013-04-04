function color_hist = computeColorHist(im,ptrs,scales,min_scale,display)
% ptrs: 2xN
% scales : 1xN

if (nargin == 4 )
    display = 0;
end

numColorHistBins = 22;
maxLValue = 100;

LBinSize = maxLValue / numColorHistBins;
abBinSize = 256 / numColorHistBins;
LBinEdges = [0:LBinSize:100];
abBinEdges = [-128:abBinSize:128];

magnif = 1.1;
scales = magnif*scales;
scales(scales < min_scale) = min_scale;

if ( ndims(im) == 2 )
    im = repmat(im,[1,1,3]);
end
% max_dim = 480;
% img_size = size(im);
% img_size = img_size(1:2);
% if ( max(img_size) > max_dim )
%     im = imresize(im, max_dim/max(img_size));
% end

[L,a,b] = RGB2Lab(double(im(:,:,1)), double(im(:,:,2)), double(im(:,:,3)));
if(max(L(:)) > 100 || min(L(:)) < 0)
    disp('error in L range');
end
if(max(a(:)) > 128 || min(a(:)) < -128)
    disp('error in a range');
end
if(max(b(:)) > 128 || min(b(:)) < -128)
    disp('error in b range\n');
end

if ( size(ptrs,2) == 2 )
    ptrs = ptrs';
end

nbins = numel(LBinEdges) + 2*numel(abBinEdges);
color_hist = zeros(nbins,size(ptrs,2));
w = size(im,2);
h = size(im,1);
for i = 1 : size(ptrs,2)
    s = scales(i);
    lx = floor( ptrs(1,i) - s + 0.5);
    rx = floor( ptrs(1,i) + s + 0.5);
    ty = floor( ptrs(2,i) - s + 0.5);
    dy = floor( ptrs(2,i) + s + 0.5);
    lx = max([lx,1]);
    rx = min([rx,w]);
    ty = max([ty,1]);
    dy = min([dy,h]);
    [x y] = meshgrid([lx:rx], [ty:dy]);
    inds = sub2ind(size(im), y(:),x(:));
    color_hist(:,i) =...
        [histc(L(inds), LBinEdges);...
         histc(a(inds), abBinEdges);...
         histc(b(inds), abBinEdges)];
    color_hist(:,i) = color_hist(:,i)./sum(color_hist(:,i));
end

if (display)
    addpath('../pwmetric/');
    while (1)
        imshow(im);
        hold on;
        plot(ptrs(1,:),ptrs(2,:),'*');
        [x y button] = ginput(1);
        if ( button == 3 )
            close;
            break;
        end
        D = slmetric_pw([x y]', ptrs, 'eucdist');
        [unused min_idx] = min(D);
        hD = slmetric_pw(color_hist(:,min_idx), color_hist, 'eucdist');
        hD(min_idx) = inf;      
        [unused min_idx2] = min(hD);

        rectangle('position',...
                [ptrs(1,min_idx)-scales(min_idx), ptrs(2,min_idx)-scales(min_idx),2*scales(min_idx), 2*scales(min_idx)],...
                'edgecolor', [1,1,0], 'linewidth', 2.0, 'curvature',[1 1]);
        rectangle('position',...
                [ptrs(1,min_idx2)-scales(min_idx2), ptrs(2,min_idx2)-scales(min_idx2),2*scales(min_idx2), 2*scales(min_idx2)],...
                'edgecolor', [1,1,0], 'linewidth', 2.0, 'curvature',[1 1]);
        hold off;
        pause;
    end    
end