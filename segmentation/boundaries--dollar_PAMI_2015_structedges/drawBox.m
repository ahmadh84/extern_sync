function [V] = drawBox(bbs, idx, segIds, edgeArr)

box = bbs(idx,:);

sz = size(segIds);
h = sz(1);  w = sz(2);

box([1 2]) = box([1 2]) - 1;

r1 = clamp(box(2) + box(4), 0, h-1); 
r0 = clamp(box(2), 0, h-1);
% box(2) = r0;
c1 = clamp(box(1) + box(3), 0, w-1); 
c0 = clamp(box(1), 0, w-1);
% box(1) = c0;

V = zeros(h, w);
validIdsmask = segIds > 0;
i = segIds(validIdsmask) + 1;

edgeso = edgeArr(:,idx);
V(validIdsmask) = edgeso(i);
end


function [v] = clamp(v, ll, ul)
if v < ll
    v = ll;
elseif v > ul;
    v = ul;
end
end
