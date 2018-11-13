function MeONSGAII(Global)
% <algorithm> <H-N>
% A Meta-Objective Approach for Many-Objective Evolutionary Optimization
% alt ---  3 --- Alternatives to the convergence component
% eta --- 1 --- Scaling Factor
% p --- 2 --- Lp Function
% C --- 0.22 --- parameter in Double-Rank Method
% Norm --- 0 --- Normlization

%--------------------------------------------------------------------------
% Copyright 2016-2018 Yiping Liu
% This is the code of MeO+NSGA-II proposed by Yiping Liu in "A Meta-
% Objective Approach for Many-Objective Evolutionary Optimization, 
% Evolutionary Computation, 2018, Early Access".
% Please contact {yiping0liu@gmail.com} if you have any problem.
% This code uses PlatEMO by Ye Tian et al.
%--------------------------------------------------------------------------

    %% Parameter setting
    [alt,eta,p,C,Norm] = Global.ParameterSet(3,1,2,0.22,1);
    fSmax = -Inf;
    fSmin = Inf;

    %% Generate random population
    Population = Global.Initialization();
    [~,FrontNo,CrowdDis,fSmax,fSmin] = MeOEnvironmentalSelection(Population,Global.N,alt,eta,p,C,fSmax,fSmin,Global.PF,Norm);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
        Offspring  = Global.Variation(Population(MatingPool));
        [Population,FrontNo,CrowdDis,fSmax,fSmin] = MeOEnvironmentalSelection([Population,Offspring],Global.N,alt,eta,p,C,fSmax,fSmin,Global.PF,Norm);
    end
end