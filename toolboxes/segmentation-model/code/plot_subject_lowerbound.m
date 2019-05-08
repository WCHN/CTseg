function plot_subject_lowerbound(lb, figname)

% ---------------------------------------------------------------------
% Get figure (create if it does not exist)
if nargin < 2
    figname = '(SPM) Subject Lower Bound';
end
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   
clf(f);

% ---------------------------------------------------------------------
% Plots

nfigs = 13;
nrows = floor(sqrt(nfigs));
ncols = ceil(nfigs/nrows);    

subplot(nrows,ncols, 1);
plot(lb.sum)
title('Lower Bound')
subplot(nrows,ncols, 2);
plot(lb.X)
title('Observations (E + KL)')
subplot(nrows,ncols, 3);
plot(lb.Z)
title('Responsibilities (KL)')
subplot(nrows,ncols, 4);
plot(lb.MU)
title('Means (KL)')
subplot(nrows,ncols, 5);
plot(lb.A)
title('Precisions (KL)')
subplot(nrows,ncols, 6);
plot(lb.lab)
title('Labels (E)')
subplot(nrows,ncols, 7);
plot(sum(lb.bf_reg,2))
title('ln(p(b))')
subplot(nrows,ncols, 8);
plot(lb.v_reg)
title('ln(p(v))')
subplot(nrows,ncols, 9);
plot(lb.aff_reg)
title('ln(p(r))')
subplot(nrows,ncols, 10);
plot(lb.lnDetbf)
title('ln|b|')
subplot(nrows,ncols, 11);
plot(lb.prop_reg)
title('ln(p(prop))')
subplot(nrows,ncols, 12);
plot(lb.mg)
title('mg')
subplot(nrows,ncols, 13);
plot(lb.ZN)
title('ZN')

drawnow
%=========================================================================