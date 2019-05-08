function out = contains(varargin)
% Backward compatibility for contains
    if exist('contains','builtin')
        out = builtin('contains', varargin{:});
        return
    end

    out = strfind(varargin{:});    
end