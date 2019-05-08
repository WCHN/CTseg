function obj = fold_cell(obj, ind)
% FORMAT obj = fold_cell(obj, ind)
% obj   - linear cell object
% ind   - Indices correspondance (obtained from unfold_cell)
%
% Transform a liner cell into a cell of cells (of cells ...).

    obj = fold(obj, ind);

end

function ind = fold(obj, ind)

    for n=1:numel(ind)
        if iscell(ind{n})
            ind{n} = fold(obj, ind{n});
        else
            ind{n} = obj{ind{n}};
        end
    end

end