function [elements, aCount] = setdiff( a, b )
a = a(:);
if isempty(a)
    elements = [];
    aCount = [];
    return
end
if isempty(b)
    elements = a;
    aCount = [1:numel(a)]';
    return
end
b = b(:);



%%

offset = min(min(a), min(b))-1;
a = a-offset;
b = b-offset;
maxVal = max(max(a), max(b));

aHist = accumarray(a, 1, [maxVal, 1]);
bHist = accumarray(b, 1, [maxVal, 1]);

elements = find(aHist>0 & bHist==0);
aCount = aHist(elements);

elements = elements+offset;


end

