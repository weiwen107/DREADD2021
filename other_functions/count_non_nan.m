function [nnan,ind] = count_non_nan(data)
% This function counts the number of non-nan elements in an array and
% returns this number and the index of non-nan elements in this array.
%   data- must be an array
%   nnan- number of non-nan elements
%   ind- the index of each non-nan elements
nnan = 0;
ind = zeros(1,numel(data));

for i = 1:numel(data)
    if isnan(data(i)) == 0
        nnan = nnan+1;
        ind(nnan) = i;
    end
end

ind = nonzeros(ind);

end

