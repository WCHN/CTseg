function [lb,gain] = check_convergence(type, lb, iter, verbose, armijo, eul_its)
% FORMAT [lb,gain] = check_convergence(type, lb, [iter], [verbose], [armijo], [eul_its])
% type    - Identifier of the current module
% lb      - Lower bound structure
% iter    - Current iteration [0]
% verbose - Verbosity level (>= 0) [0]
% armijo  - Armijo factor for which success was reached (gauss-newton only)
% eul_its - Number of integration steps used (velocities only)
%
% Compute lower bound (by summing its parts) and its gain
% + print info
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 3, iter    = 0; end
if nargin < 4, verbose = 0; end
if nargin < 5, armijo  = []; end
if nargin < 6, eul_its = []; end

lb.sum(end + 1) =   lb.X(end) + lb.lnDetbf(end) + lb.Z(end) ...
                  + lb.MU(end) + lb.A(end) + sum(lb.bf_reg(end,:)) ...
                  + lb.v_reg(end) + lb.lab(end) + lb.aff_reg(end) ...
                  + lb.prop_reg(end) + lb.mg(end) + lb.ZN(end);

if verbose >= 2
    plot_subject_lowerbound(lb)
end

gain = spm_misc('get_gain',lb.sum);

if verbose >= 1
    if     numel(lb.sum) < 2  
        incr = '';
    elseif lb.sum(end) > lb.sum(end-1) 
        incr = '(+)';
    elseif lb.sum(end) < lb.sum(end-1) && gain > 1e-6
        incr = '(-)';
    else                         
        incr = '(=)';
    end
    fprintf('%3s | %3d | lb = %10.6f | gain = %1.5f | %3s', type, iter, lb.sum(end), gain, incr);
    
    if ~isempty(armijo)
        fprintf(' | armijo = ');
        for c=1:numel(armijo)
            fprintf('%1.5f ', armijo(c));
        end
    end
    
    if ~isempty(eul_its)
        fprintf(' | eul_its = %i ', eul_its);        
    end

    fprintf('\n');
end
lb.last = lb.sum(end);
% =========================================================================