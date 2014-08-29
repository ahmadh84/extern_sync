function [ elements, aCount, bCount ] = intersect( a, b )

if isempty(a) || isempty(b)
    elements = [];
    aCount = [];
    bCount = [];
    return
end

a = a(:);
b = b(:);

%%
offset = min(min(a), min(b))-1;
a = a-offset;
b = b-offset;
maxVal = max(max(a), max(b));

aHist = accumarray(a, 1, [maxVal, 1]);
bHist = accumarray(b, 1, [maxVal, 1]);

elements = find(aHist>0 & bHist>0);
aCount = aHist(elements);
bCount = bHist(elements);

elements = elements+offset;

end

