function B = get_rigid_basis(is3d)
% Basis functions for the lie algebra of the special Eucliden group
% (SE(3)): translation and rotation.
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 1, is3d = false; end

if is3d
    B                = zeros(4,4,6);
    B(1,4,1)         = 1;
    B(2,4,2)         = 1;
    B(3,4,3)         = 1;
    B([1,2],[1,2],4) = [0 1;-1 0];
    B([3,1],[3,1],5) = [0 1;-1 0];
    B([2,3],[2,3],6) = [0 1;-1 0];
else
    B                = zeros(4,4,3);
    B(1,4,1)         = 1;
    B(2,4,2)         = 1;
    B([1,2],[1,2],3) = [0 1;-1 0]; 
end
%==========================================================================