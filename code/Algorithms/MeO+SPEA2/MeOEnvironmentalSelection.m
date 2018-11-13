function [Population,Fitness,fSmax,fSmin] = MeOEnvironmentalSelection(Population,N,alt,eta,p,C,fSmax,fSmin,TPF,Norm)
% The MeO environmental selection of SPEA2
% The code is modified from EnvironmentalSelection.m of SPEA2 in PlatEMO

%--------------------------------------------------------------------------
% Copyright 2016-2018 Yiping Liu
% This is the code of MeO+SEPA2 proposed by Yiping Liu in "A Meta-
% Objective Approach for Many-Objective Evolutionary Optimization, 
% Evolutionary Computation, 2018, Early Access".
% Please contact {yiping0liu@gmail.com} if you have any problem.
% This code uses PlatEMO by Ye Tian et al.
%--------------------------------------------------------------------------

    %% Meta Objective
    [MetaObj,fSmax,fSmin] = MetaObjective(Population.objs,alt,eta,p,C,fSmax,fSmin,TPF,Norm);

    %% Calculate the fitness of each solution
    Fitness = CalFitness(MetaObj);

    %% Environmental selection
    Next = Fitness < 1;
    if sum(Next) < N
        [~,Rank] = sort(Fitness);
        Next(Rank(1:N)) = true;
    elseif sum(Next) > N
        Del  = Truncation(MetaObj(Next,:),sum(Next)-N);
        Temp = find(Next);
        Next(Temp(Del)) = false;
    end
    % Population for next generation
    Population = Population(Next);
    Fitness    = Fitness(Next);
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    %% Truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end