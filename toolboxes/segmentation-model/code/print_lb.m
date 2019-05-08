function print_lb(lb,it)
if nargin < 2, it = 0; end

fields = fieldnames(lb);
for i=1:numel(fields)
    print_str(lb,fields{i},it);
end
fprintf('\n');
%==========================================================================

%==========================================================================
function print_str(lb,field,it)
val = lb.(field);
if numel(val) > it + 1
    fprintf('| %s = %10.6f ',field,val(end - it));
else
    fprintf('| %s = %10.6f ',field,val(end));
end
%==========================================================================