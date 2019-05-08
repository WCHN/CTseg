function Z = cluster2template(Z_cluster,part)
% FORMAT Z = cluster2template(Z_cluster,part)
% Z_cluster - GMM cluster responsibilities
% part      - Stucture with lkp field that maps clusters to classes.
% Z         - Class responsibilities
%
% Create class-responsibilities from GMM cluster responsibilities.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
Iz  = size(Z_cluster,1);
lkp = part.lkp;
K   = max(lkp);
Z   = zeros([Iz K]);
for k=1:K
    Z(:,k) = sum(Z_cluster(:,lkp == k),2);
end
%==========================================================================