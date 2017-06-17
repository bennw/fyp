%%%%%%%%% PARAMETER CONTROL %%%%%%%%%
HighMovementEntropy = 0;        % 1 for high entropy, 0 for low entropy
gamma = 0.80;                   % discount factor
PlotLenVsCost = 0;              % 1 to plot mean cost vs. path length
%%%%%%% END PARAMETER CONTROL %%%%%%%

% TODO: explore relationship between gamma and improvement% with pathlen

BSize = 6;              % b = 0..BSize-1
N = 9;                  % parameter N
joules = 8.85;
c_fix = joules;         % parameter gamma^b
c_c = [3 5 3.5 4 6 4 3.5 5 0]*joules; % parameter gamma^pc(i)
C = zeros(N,5);         % cost matrix C
P = [40 30 0 0 30;
     40 20 0 20 20;
     40 0 0 30 30;
     40 20 20 0 20;
     40 15 15 15 15;
     40 0 20 20 20;
     40 30 30 0 0;
     40 20 20 20 0;
     100 0 0 0 0]/100;
if HighMovementEntropy == 0
    P = [10 80 0 0 10;
         6 80 0 6 8;
         10 0 0 10 80;
         6 8 80 0 6;
         80 5 5 5 5;
         6 0 6 80 8;
         10 10 80 0 0;
         6 6 8 80 0;
         100 0 0 0 0]/100;
end
V = zeros(BSize*N, 5);
Vold = zeros(BSize*N, 5);
Dd = zeros(BSize*N, 5); % decisionparameter D^d
DdOld = zeros(BSize*N, 5);
alphaintermediate = 0;
alpha = [-10000 alphaintermediate 10000];
lambda = 2;             % parameter lambda
epsilon = 0.1*c_fix;    % parameter epsilon (in termination condition)
trials = 100000; % number of trials in test phase
pathLen = 60;           % path length in test phase
% cost: matrix recording mean costs obtained during trials for plotting
%       graphs. Size of first dimension: number of caching setups (e.g.
%       size 3 with reactive, optimal and pure proactive algorithms). Size
%       of second dimension: 3 - cost(:,2) contians mean caching costs,
%       cost(:,3) contains mean download costs.
cost = zeros(size(alpha,2),3);
cost(:,1) = alpha.';
if PlotLenVsCost == 1
    lenVsCost = zeros(size(alpha,2),pathLen);
    lenVsCostFreq = zeros(size(alpha,2),pathLen);
end

% 1) Populate trial paths 
ssDist = [1,0,0,0,0,0,0,0,0];
cell = zeros(trials,pathLen+1);
dirarr = zeros(trials,pathLen);
for path=1:trials
    [cell(path,:), dirarr(path,:)] = generatePath(P, ssDist, pathLen);
end

for alph=1:size(alpha,2)
    % 2) Compute matrix C and caching decision matrix Dc
    for i=1:9
        for dir=0:4
            j = movecell(i,dir);
            if j == -1 || j == 9
                C(i,dir+1) = 0;
                Dc(i, dir+1) = 0;
            else
                c_j = c_c(j);
                a = P(i,dir+1)*c_j + c_fix;         % expected cost of caching
                b = P(i,dir+1)*(lambda*c_j + c_fix); % expected cost of no cache
                if a < alpha(alph)+b && i ~= 9
                    C(i,dir+1) = c_j;
                    Dc(i, dir+1) = 1;
                else
                    C(i,dir+1) = lambda*c_j + c_fix;
                    Dc(i, dir+1) = 0;
                end
            end
        end
    end

    Cp = C.';

    % 3) infinite horizon iteration
    V = zeros(BSize*N, 5);
    for iter=1:2000
        cp  = repmat( sum(Dc(:,:),2)*c_fix, 1, 5);
        cpp = C;
        Vold = V;
        Ddold = Dd;
        pvSum = sum(repmat(P,BSize,1) .* Vold, 2); % vector with B*9 (B*n) cells
        for b=0:BSize-1
            if b == 0
                V(N*b+1:N*b+N,:) = gamma*jVec_to_iDirArr( pvSum(N*(b+1)+1:N*(b+1)+N) ) + cp + cpp;
                Dd(N*b+1:N*b+N,:) = 1;
            elseif b >= BSize-1
                V(N*b+1:N*b+N,:) = gamma*jVec_to_iDirArr( pvSum(N*(b-1)+1:N*(b-1)+N) ) + cp;
                Dd(N*b+1:N*b+N,:) = 0;
            else
                costZero = gamma*jVec_to_iDirArr( pvSum(N*(b-1)+1:N*(b-1)+N) ) + cp;
                costOne = gamma*jVec_to_iDirArr( pvSum(N*(b+1)+1:N*(b+1)+N) ) + cp + cpp;
                V(N*b+1:N*b+N,:) = min(costZero, costOne);
                Dd(N*b+1:N*b+N,:) = V(N*b+1:N*b+N,:) ~= costZero;
            end
        end
        residue = max(max(V - Vold).');
        if residue < epsilon*(1-gamma)/(2*gamma)
            break;
        end
    end
    
    % 4) Feed Dd back into new Dc (Dc2)
    % Dc2 = 0 if Dd = 0
    Dc2 = repmat(Dc,BSize,1) .* Dd;

    % 5) Test with trial paths from part 1
    for path=1:trials
        buf = 0; % current user buffer
        for ind=1:pathLen
            % 5.0: path ends if user reaches home cell
            if cell(path,ind) == 9
                break
            end
            % 5.1: find out if user chooses to download
            downloadDecision = Dd(buf*N + cell(path,ind),dirarr(path,ind));
            if buf == 0 % no user buffer left, must download
                downloadDecision = 1;
            elseif buf >= BSize-1 % no need to download anymore
                downloadDecision = 0;
            end
            % 5.2: compute costs
            OCP = Dc2( N*buf + 1 : N*buf + N,:); % optimal caching policy - 9x5 matrix corresponding to i,dir(i,j)
            if OCP(cell(path,ind),dirarr(path,ind)) == 0
                multiplier = lambda;
                if downloadDecision == 1
                    cost(alph,2) = cost(alph,2) + c_fix;
                end
            else
                multiplier = 1;
            end
            cost(alph,2) = cost(alph,2) + sum(OCP(cell(path,ind),:),2)*c_fix;  % caching cost (2)
            cost(alph,3) = cost(alph,3) + downloadDecision*multiplier*c_c(cell(path,ind+1)); % download cost (3)
            if downloadDecision == 1
                buf = buf + 1;
            else
                buf = buf - 1;
            end
        end
        % START OF LEN VS COST MODIFICATION
        if PlotLenVsCost == 1
            if ind <= 10
                lenVsCost(alph,1) = lenVsCost(alph,1) + cost(alph,3);
                lenVsCostFreq(alph,1) = lenVsCostFreq(alph,1) + 1;
            elseif ind <= 20
                lenVsCost(alph,2) = lenVsCost(alph,2) + cost(alph,3);
                lenVsCostFreq(alph,2) = lenVsCostFreq(alph,2) + 1;
            elseif ind <= 30
                lenVsCost(alph,3) = lenVsCost(alph,3) + cost(alph,3);
                lenVsCostFreq(alph,3) = lenVsCostFreq(alph,3) + 1;
            elseif ind <= 40
                lenVsCost(alph,4) = lenVsCost(alph,4) + cost(alph,3);
                lenVsCostFreq(alph,4) = lenVsCostFreq(alph,4) + 1;
            elseif ind <= 50
                lenVsCost(alph,5) = lenVsCost(alph,5) + cost(alph,3);
                lenVsCostFreq(alph,5) = lenVsCostFreq(alph,5) + 1;
            end
            cost(alph,2) = 0;
            cost(alph,3) = 0;
        end
        % END OF LEN VS COST MODIFICATION
    end
end
cost(:,2:3) = cost(:,2:3) / trials;

% 6) Plots


if PlotLenVsCost == 1
    hold on;
    lenVsCost = lenVsCost./lenVsCostFreq;
    plot(lenVsCost(2,1:5),'^-');
    legend('Gamma = 0.5', 'Gamma = 0.7', 'Gamma = 0.9')
    xlabels = {'4-9'; '10-19'; '20-29'; '30-39'; '40-49' };
    set(gca,'xtick',1:5)
	set(gca,'xticklabel',xlabels)
    xlabel('Path length');
    ylabel('Mean cost');
    if HighMovementEntropy == 0
        title('Infinite horizon, low movement entropy');
    else
        title('Infinite horizon, high movement entropy');
    end
else
    plot(cost(:,2),'^-');
    hold on
    plot(cost(:,3),'d-');
    plot(cost(:,2)+cost(:,3),'*-');
    legend('Caching cost','Download cost','Caching + download cost')
    xlabels = {'reactive'; 'optimal'; 'pure proactive' };
    set(gca,'xtick',[1 2 3])
    set(gca,'xticklabel',xlabels)
    xlabel('Caching algorithm');
    ylabel('Mean cost');
    if HighMovementEntropy == 0
        title('Infinite horizon, low movement entropy');
    else
        title('Infinite horizon, high movement entropy');
    end
end