%==========================================================================
function X = set_slices(X, dim, ind, S)
% FORMAT X = set_slices(X, dim, ind, S)
% X   - input array
% dim - dimension(s) along which to select a slice
% ind - indices to select in each slice
% S   - slice
sizeX    = size(X);
sub      = struct;
sub.type = '()';
sub.subs = repmat({':'}, [1 numel(sizeX)]);
if iscell(ind)
    assert(numel(dim) == numel(ind));
    for d=1:numel(dim)
        sub.subs{dim(d)} = ind{d};
    end
else
    sub.subs{dim} = ind;
end
X = subsasgn(X,sub,S);
%==========================================================================
