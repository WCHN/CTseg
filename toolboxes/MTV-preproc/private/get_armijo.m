function armijo = get_armijo(dat)
% Rigid optimisation line-search parameter
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

C      = numel(dat);
armijo = cell(1,C);
for c=1:C   
    armijo{c} = ones(1,dat(c).N);
end
%==========================================================================