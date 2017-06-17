function [ cell, dirarr ] = generatePath( P, startPosDist, pathLen )
%GENERATEPATH Summary of this function goes here
%   P: cell transition matrix
%   ssDist: distribution function (probability) of starting cell
%   pathlen: number of cell transitions (cells visited = pathLen+1)
%   OUTPUT: dirarr: range [1,5]

    i = randarr(startPosDist);
    cell = zeros(1,pathLen+1);
    dirarr = zeros(1,pathLen); % possible values in dirarr: [1,5]
    cell(1) = i;
    % propagate path
    for iter=1:pathLen
        dir = randarr(P(cell(iter),:));
        dirarr(iter) = dir;
        cell(iter+1) = movecell(cell(iter), dir-1);
    end
end

