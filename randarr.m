function [ ind ] = randarr( arr )
% for stochastic array arr, return index with probability arr(ind)

ind = -1;
r = rand;
for i=1:numel(arr)
    r = r - arr(i);
    if r < 0
        ind = i;
        break
    end
    ind = i;
end

end

