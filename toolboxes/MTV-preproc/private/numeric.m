function a = numeric(a)
% FORMAT a = numeric(fa)
% fa - file_array (or numeric array)
% a  - numeric array
% 
% Convert to numeric form.
%   > Wrapper around 'file_array>numeric' so that it also works with 
%     numeric arrays (in this case it does nothing).

    if isa(a, 'file_array')
        a = numeric(a);
    end
    
end