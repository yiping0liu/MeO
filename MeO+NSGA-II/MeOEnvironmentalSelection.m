function [Population,FrontNo,CrowdDis,fSmax,fSmin] = MeOEnvironmentalSelection(Population,N,alt,eta,p,C,fSmax,fSmin,TPF,Norm)
% The MeO environmental selection of NSGA-II
% The code is modified from EnvironmentalSelection.m of NSGA-II in PlatEMO

%--------------------------------------------------------------------------
% Copyright Yiping Liu
% This is the code of MeO+NSGA-II proposed by Yiping Liu in "A Meta-
% Objective Approach for Many-Objective Evolutionary Optimization, 
% Evolutionary Computation, 2020, 28(1), pp.1-25".
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

    %% Meta Objective
    [MetaObj,fSmax,fSmin] = MetaObjective(Population.objs,alt,eta,p,C,fSmax,fSmin,TPF,Norm);

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(MetaObj,N);
    Next = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(MetaObj,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end