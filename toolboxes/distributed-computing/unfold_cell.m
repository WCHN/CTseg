function [obj, ind] = unfold_cell(obj, depth)
% FORMAT [obj, ind] = unfold_cell(obj, (depth))
% obj   - cell of cells object
% depth - Number of levels to unfold [2]
% ind   - Indices correspondance
%
% Transform a cell of cells (of cells ...) into a linear cell.

    if nargin < 2
        depth = 2;
    end

    [obj, ind] = unfold(obj, 1, depth);

end

% function N = all_numel(obj, d, maxd)
%     N = 0;
%     if d < maxd
%        for n=1:numel(obj)
%            N = N + all_numel(obj{n}, d+1, maxd);
%        end
%     elseif d == maxd
%         N = numel(obj);
%     end
% end

function [objout, ind, curn] = unfold(obj, d, maxd, curn)

    if nargin < 4
        curn = 1;
    end

    ind = cell(size(obj));
    N1  = numel(obj);
    
    if d < maxd
        objout = {};
        for n=1:N1
            ndim = numel(size(obj));
            [pos{1:ndim}] = ind2sub(size(obj), n);
            [obj1, ind1, curn] = unfold(obj{n}, d+1, maxd, curn);
            objout = [objout obj1];
            ind = subsasgn(ind, struct('type', '{}', 'subs', {pos}), ind1);
        end
    elseif d == maxd
        objout = obj;
        ndim   = numel(size(obj));
        for n=1:N1
            [pos{1:ndim}] = ind2sub(size(obj), n);
            ind = subsasgn(ind, struct('type', '{}', 'subs', {pos}), curn);
            curn = curn + 1;
        end
        1;
    else
        error('Should not have reached this depth')
    end
    
end