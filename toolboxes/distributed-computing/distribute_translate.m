function varargout = distribute_translate(opt, varargin)
% FORMAT [out1, ...] = distribute_translate(opt, arg1, ...,)
%
% opt - Option structure with fields:
%     translate - Cell array of size 2xN with translation between client 
%                 and server paths [{client.folder server.folder}].
%                 Example:
%                     {'/home/me/'     '/mnt/users/me/' ;
%                      '/shared/data/' '/mnt/shared/data'}
%     restrict  - Restrict translation to either 'char' or 'file_array'
%                 objects. If empty, do not restrict.
% arg   - Argument that needs translation (string, file_array, struct, ...)

    if isempty(opt.translate)
        varargout = varargin;
    end

    if isfield(opt, 'server_to_client') && opt.server_to_client
        if ~isempty(opt.translate)
            opt.translate = {opt.translate{:,2} opt.translate{:,1}};
        end
    else
        opt.server_to_client = false;
    end

    for i=1:min(numel(varargin), nargout)
        varargout{i} = translate(opt, varargin{i});
    end
    
end

function obj = translate(opt, obj)
    
    if ischar(obj) && ~strcmpi(opt.restrict, 'file_array')
        [pth,nam,ext] = fileparts(obj);
        if strcmp(ext,'.mat')
            if opt.server_to_client
                pth = translate(opt, fullfile(pth,nam));
                obj = [char(pth) ext];
            end
            obj1 = load(obj);
            obj1 = translate(opt, obj1);
            save(obj,'-struct','obj1');
            clear obj1
            if ~opt.server_to_client
                pth = translate(opt, fullfile(pth,nam));
                obj = [char(pth) ext];
            end
        else
            for j=1:size(opt.translate, 1)
                obj = strrep(obj, opt.translate{j,1}, opt.translate{j,2});
            end
        end
        return
    end
    if isa(obj, 'file_array') && ~strcmpi(opt.restrict, 'char')
        obj.fname = translate(opt, obj.fname);
        return
    end
    if isa(obj, 'nifti')
        N = numel(obj);
        for n=1:N
            obj(n).dat = translate(opt, obj(n).dat);
        end
        return
    end
    if isstruct(obj)
        fields = fieldnames(obj);
        for j=1:numel(fields)
            field = fields{j};
            for n=1:numel(obj)
                obj(n).(field) = translate(opt, obj(n).(field));
            end
        end
        return
    end
    if iscell(obj)
        for j=1:numel(obj)
            obj{j} = translate(opt, obj{j});
        end
        return
    end
    if isa(obj, 'matlab.ui.Figure')
        obj = [];
        return
    end

end