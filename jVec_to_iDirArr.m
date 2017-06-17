function [ iDirArr ] = jVec_to_iDirArr( v )
% Populates a 9x5 matrix (i x dir(i,j)) matrix from input vector v.
% The value from the jth index of v will be populated into all cells in
% iDirArr with dir(i,j) corresponding to j
%   input: 9x1 vector v denoting cell j
%   output: 9x5 array iDirArr denoting i x j, where i runs from 1-9 in
%   first column

iDirArr = zeros(9, 5);
for i=1:9
    for j=1:5
        jCell = movecell(i,j-1);
        if jCell == -1
            iDirArr(i,j) = 0;
        else
            iDirArr(i,j) = v(jCell);
        end
    end
end