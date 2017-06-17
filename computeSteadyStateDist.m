function [ dist ] = computeSteadyStateDist( P )
%COMPUTESTEADYSTATEDIST Summary of this function goes here
%   Compute steady state distribution of a 9 cell space given transition
%   probabilities P(i, dir(i,j))
i = 1;
n = 1000000; % number of iterations
freq = [0 0 0 0 0 0 0 0 0];
for iter=1:n
    dir = randarr(P(i,:));
    i = movecell(i, dir-1);
    freq(i) = freq(i) + 1;
end
dist = freq/n;

end

