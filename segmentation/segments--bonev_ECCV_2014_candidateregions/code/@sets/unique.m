function [elements, counts] = unique(a)

a = a(:);

offsetA = min(a)-1;
a = a-offsetA;
counts = accumarray(a, 1);
elements = find(counts);
counts = counts(elements);
elements = elements+offsetA;

end

