function [ elements, aCount, bCount ] = union( a, b )
a = a(:);
b = b(:);

if isempty(a)
    if isempty(b)
        elements = [];
        aCount = [];
        bCount = [];
    else
        elements = b;
        aCount = [];
        bCount = [1:numel(b)]';
    end
    return
end
if isempty(b)
    elements = a;
    aCount = [1:numel(a)]';
    bCount = [];
    return
end


%%

offset = min(min(a), min(b))-1;
a = a-offset;
b = b-offset;
maxVal = max(max(a), max(b));

aHist = accumarray(a, 1, [maxVal, 1]);
bHist = accumarray(b, 1, [maxVal, 1]);

elements = find(aHist>0 | bHist>0);
aCount = aHist(elements);
bCount = bHist(elements);

elements = elements+offset;


end

