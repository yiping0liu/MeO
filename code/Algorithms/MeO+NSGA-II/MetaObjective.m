function [MetaObj,fSmax,fSmin] = MetaObjective(OriObj,alt,eta,p,C,fSmax,fSmin,TPF,Norm)
% Meta Objective Calculation

%--------------------------------------------------------------------------
% Copyright 2016-2018 Yiping Liu
% This is the code of MeO+NSGA-II proposed by Yiping Liu in "A Meta-
% Objective Approach for Many-Objective Evolutionary Optimization, 
% Evolutionary Computation, 2018, Early Access".
% Please contact {yiping0liu@gmail.com} if you have any problem.
% This code uses PlatEMO by Ye Tian et al.
%--------------------------------------------------------------------------

[N,M] = size(OriObj);

%% Normalization
if alt<4
   Zmin = min(OriObj,[],1);
   if Norm
       Zmax = max(OriObj,[],1);
       OriObj = (OriObj-repmat(Zmin,N,1))./(repmat(Zmax,N,1)-repmat(Zmin,N,1));
   else
       OriObj = OriObj-repmat(Zmin,N,1);
   end
else
   Zmin = min(TPF,[],1);
   if Norm
       Zmax = max(TPF,[],1);
       OriObj = (OriObj-repmat(Zmin,N,1))./(repmat(Zmax,N,1)-repmat(Zmin,N,1));
       TPF = (TPF-repmat(Zmin,size(TPF,1),1))./(repmat(Zmax,size(TPF,1),1)-repmat(Zmin,size(TPF,1),1));
   else
       OriObj = OriObj-repmat(Zmin,N,1);
       TPF = TPF-repmat(Zmin,size(TPF,1),1);
   end
end

%% Convergence Component
switch alt    
    case 1 %Pareto Rank
        [Rank,~] = NDSort(OriObj,Inf);
        fC = Rank'-1; 
        
    case 2 %L_p Function
        fS = sum(OriObj.^p,2).^(1/p); 
        fSmin = min(fSmin,min(fS));
        fSmax = max(fSmax,max(fS));
        fC = (fS-fSmin)./(fSmax-fSmin); 
        
    case 3 %Double Rank
        [Rank,~] = NDSort(OriObj,Inf);
        fS = sum(OriObj.^p,2).^(1/p);
        nMax = 0;
        Rankl = zeros(N,1);
        for i=1:N
            Cosine = 1 - pdist2(OriObj(i,:),OriObj,'cosine');
            Angle = acos(Cosine);
            Neighbours = find(Angle<C*pi*N^(-1/(M-1)));
            nMax = max(nMax,length(Neighbours));
            [~,I] = sort(fS(Neighbours));
            Rankl(i) = find(I==find(Neighbours==i));
        end
        fC = (Rank'-1).*2+(Rankl-1)./nMax;
        
    case 4 %Distance to True PF
        fC = min(pdist2(OriObj,TPF),[],2);

end

%% Diversity Component
fD = (1 - pdist2(OriObj,eye(M),'cosine')).^2;

%% Meta Objective
MetaObj = fD+repmat(fC,1,M).*eta;


end