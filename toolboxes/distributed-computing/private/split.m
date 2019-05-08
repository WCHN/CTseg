function out = split(varargin)
% Backward compatibility for split
    if exist('split','builtin')
        out = builtin('split', varargin{:});
        return
    end

    out = strsplit(varargin{:});    
end