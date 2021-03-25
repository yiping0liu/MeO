function [Population,Fitness,fSmax,fSmin] = MeOEnvironmentalSelection(Population,N,alt,eta,p,C,fSmax,fSmin,TPF,Norm)
% The MeO environmental selection of SPEA2
% The code is modified from EnvironmentalSelection.m of SPEA2 in PlatEMO

%--------------------------------------------------------------------------
% Copyright Yiping Liu
% This is the code of MeO+SPEA2 proposed by Yiping Liu in "A Meta-
% Objective Approach for Many-Objective Evolutionary Optimization, 
% Evolutionary Computation, 2020, 28(1), pp.1-25".
% Please contact {yiping0liu@gmail.com} if you have any problem.
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