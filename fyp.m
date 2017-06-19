%%%%%%%%% PARAMETER CONTROL %%%%%%%%%
PolarisedProbabilities = 1;     % 1 for polarised probabilities, 0 otherwise
HighBuffer = 0;                 % 1 for high buffer, 0 for low buffer
%%%%%%% END PARAMETER CONTROL %%%%%%%

BSize = 6; % b = 0..BSize-1
if HighBuffer == 1
    BSize = 21;
end
LSize = 61; % l = 0Size..L-1
B = BSize-1;            % parameter B
L = LSize-1;            % parameter L
N = 9;                  % parameter N
joules = 8.85;
% alphaintermediate = -.6:0.2:.4*joules;
alphaintermediate = 0;
alpha = [-10000 alphaintermediate 10000];
c_fix = joules;         % parameter gamma^b
c_c = [3 5 3.5 4 6 4 3.5 5 3]*joules; % parameter gamma^pc(i)
C = zeros(N,5);         % cost matrix C
P = [40 30 0 0 30;
     40 20 0 20 20;
     40 0 0 30 30;
     40 20 20 0 20;
     40 15 15 15 15;
     40 0 20 20 20;
     40 30 30 0 0;
     40 20 20 20 0;
     40 0 30 30 0]/100; % parameter p(i,j)
if PolarisedProbabilities == 1
    P = [10 80 0 0 10;
         6 80 0 6 8;
         10 0 0 10 80;
         6 8 80 0 6;
         80 5 5 5 5;
         6 0 6 8 80;
         10 10 80 0 0;
         6 6 8 80 0;
         10 0 10 80 0]/100;
end
Vs = zeros(N, 5*BSize*LSize);
Dd = zeros(N, 5*BSize*LSize);   % decision parameter D^d
PVs = zeros(N, BSize*LSize);
lambda = 2;                     % parameter lambda
trials = 10000;                 % number of trials in test phase
% cost: matrix recording mean costs obtained during trials for plotting
%       graphs. Size of first dimension: number of caching setups (e.g.
%       size 3 with reactive, optimal and pure proactive algorithms). Size
%       of second dimension: 3 - cost(:,2) contians mean caching costs,
%       cost(:,3) contains mean download costs.
cost = zeros(size(alpha,2),3);
cost(:,1) = alpha.';

% Compute steady state distribution
ssDist = computeSteadyStateDist(P);

for alph=1:size(alpha,2)
    % 1) Compute matrix C and caching decision matrix Dc
    for i=1:9
        for dir=0:4
            j = movecell(i,dir);
            if j == -1
                C(i,dir+1) = 0;
                Dc(i, dir+1) = 0;
            else
                c_j = c_c(j);
                a = P(i,dir+1)*c_j + c_fix;         % expected cost of caching
                b = P(i,dir+1)*(lambda*c_j + c_fix); % expected cost of no cache
                if a < alpha(alph)+b
                    C(i,dir+1) = c_j;
                    Dc(i, dir+1) = 1;
                else
                    C(i,dir+1) = lambda*c_j + c_fix;
                    Dc(i, dir+1) = 0;
                end
            end
        end
    end

    % 2) populate PX matrix (size bl x j) and dl decision matrix Dd
    % Convention: (b,l) runs from (0,0) (0,1) .. (0,L) (1,0) (1,1) ..
    cp  = repmat( sum(Dc(:,:),2)*c_fix, 1, 5);
    cpp = C;
    for l=1:L
        for b=0:B
            PVs_ind = b*LSize + l + 1;
            Vs_ind = 5*(b*LSize + l) + 1 : 5*(b*LSize + l) + 5;
            if b == 0
                Vs(:,Vs_ind) = cp + cpp + jVec_to_iDirArr( PVs(:,PVs_ind-1) );
                Dd(:,Vs_ind) = 1;
            else
                Vs(:,Vs_ind) = cp + cpp + jVec_to_iDirArr( PVs(:,PVs_ind-1) ); % Dd = 1
                DZero = cp + jVec_to_iDirArr( PVs(:,PVs_ind-LSize) ) ;             % Dd = 0
                Vs(:,Vs_ind) = min(Vs(:,Vs_ind), DZero); % V(s) = c' + min( Dd=1, Dd=0 )
                Dd(:,Vs_ind) = Vs(:,Vs_ind) ~= DZero;
            end
            PVs(:,PVs_ind) = sum(P .* Vs(:,Vs_ind), 2);
        end
    end

    % 3) Feed Dd back into new Dc (Dc2)
    
    % Dc2 = 0 if Dd = 0
    Dc2 = repmat(Dc,1,BSize*LSize) .* Dd;
    Dmat = Dd(:,5*(B*LSize + L)+1:5*(B*LSize + L)+5).';

    % 4) Make path where initial cell is randomised by SteadyStateArr
    % Compute path cost and store in AlphaCcostDcost
    for path=1:trials
        % generate path for this iteration
        [cell, dirarr] = generatePath(P, ssDist, B+L);
        % initialise temporary counters for path tracing
        b = B;
        l = L;
        for ind=1:B+L
            % 0: find out if user chooses to download
            downloadDecision = Dd(cell(ind),5*(b*LSize + l)+dirarr(ind));
            if l == 0 % no need to download anymore
                break;
            elseif b == 0 % no user buffer left, must download
                downloadDecision = 1;
            end
            OCP = Dc2( :,5*(b*LSize + l) + 1 : 5*(b*LSize + l) + 5 ).'; % optimal caching policy - 5x9 matrix corresponding to i,dir(i,j)
            if OCP(dirarr(ind),cell(ind)) == 0
                multiplier = lambda;
                if downloadDecision == 1
                    cost(alph,2) = cost(alph,2) + c_fix; % caching cost incurred during reactive download
                end
            else
                multiplier = 1;
            end
            cost(alph,3) = cost(alph,3) + downloadDecision*multiplier*c_c(cell(ind+1));  % download cost
            cost(alph,2) = cost(alph,2) + sum(OCP(:,cell(ind)),1)*c_fix;                 % caching cost
            if downloadDecision == 1
                l = l - 1;
            else
                b = b - 1;
            end
        end
    end
end

cost(:,2:3) = cost(:,2:3) / trials;
polycoeff = polyfit(cost(:,2), cost(:,3), 4);
polyx = 20:1:72;
polyy = polyval(polycoeff, polyx);

% f1 = figure;
% 
% bar(disp(:,:), 'stacked')
% barxlabels = alpha;
% set(gca,'xticklabel',barxlabels)
% ylabel('Expected cost');
% title('Expected cost vs. alpha, B=3, L=18');

% f2 = figure;
% 
% bar(alphaCcostDcost(:,[3 2]), 'stacked')
% barxlabels = {'reactive'; 'optimal'; 'pure proactive' };
% set(gca,'xticklabel',barxlabels)
% ylabel('Mean cost');
% title('Mean cost vs. alpha, B=3, L=18, 10000 simulations with random paths');
% 
% f3 = figure;
% %x = [mean(reactiveCost(1,:),2), mean(ourAlgoCost(1,:),2), mean(proactiveCost(1,:),2)];
% %y = [mean(reactiveCost(2,:),2), mean(ourAlgoCost(2,:),2), mean(proactiveCost(2,:),2)];
% hold on
% plot(alphaCcostDcost(:,2), alphaCcostDcost(:,3), 'x-');
% plot(polyx, polyy, 'r--')
% hold off
% xlabel('Mean caching cost');
% ylabel('Mean download cost');
% title('Download vs. caching cost comparison, B=3, L=18, 10000 simulations with random paths');
% for i=1:size(alpha,2)
%     text(alphaCcostDcost(i,2)+0.25, alphaCcostDcost(i,3)+0.4, num2str(alpha(i)));
% end

hold on
%plot(cost(:,2),'*-');
%plot(cost(:,3),'*-');
plot(cost(:,2)+cost(:,3),'*-');
legend('Scenario 1.1: low buffer','Scenario 1.2: high buffer','Scenario 1.3: low buffer, polarised probabilities')
xlabels = {'reactive'; 'optimal'; 'pure proactive' };
set(gca,'xtick',[1 2 3])
set(gca,'xticklabel',xlabels)
xlabel('Caching algorithm');
ylabel('Sum of mean caching and download costs');
title('Finite horizon, total cost');
