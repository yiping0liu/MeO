classdef MeOSPEA2 < ALGORITHM
% <multi/many> <real/binary/permutation>
% A Meta-Objective Approach for Many-Objective Evolutionary Optimization
% alt ---  3 --- Alternatives to the convergence component
% eta --- 1 --- Scaling Factor
% p --- 2 --- Lp Function
% C --- 0.22 --- parameter in Double-Rank Method
% Norm --- 0 --- Normlization

%--------------------------------------------------------------------------
% Copyright Yiping Liu
% This is the code of MeO+SPEA2 proposed by Yiping Liu in "A Meta-
% Objective Approach for Many-Objective Evolutionary Optimization, 
% Evolutionary Computation, 2020, 28(1), pp.1-25".
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [alt,eta,p,C,Norm] = Algorithm.ParameterSet(3,1,2,0.22,1);
            fSmax = -Inf;
            fSmin = Inf;

            %% Generate random population
            Population = Problem.Initialization();
            [MetaObj,fSmax,fSmin] = MetaObjective(Population.objs,alt,eta,p,C,fSmax,fSmin,Problem.optimum,Norm);
            Fitness    = CalFitness(MetaObj);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,Fitness);
                Offspring  = OperatorGA(Population(MatingPool));
                [Population,Fitness,fSmax,fSmin] = MeOEnvironmentalSelection([Population,Offspring],Problem.N,alt,eta,p,C,fSmax,fSmin,Problem.optimum,Norm);
            end
        end
    end
end