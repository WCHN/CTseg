function varargout = distribute(varargin)
% _________________________________________________________________________
%
%             Distribute Matlab jobs locally or on a cluster
%
% -------------------------------------------------------------------------
% ARGUMENTS
% ---------
%
% FORMAT [opt,out1, ...] = distribute(opt, func, ('iter'/'inplace'), arg1, ...)
%
% opt  - Option structure. See 'help distribute_default'.
% func - Matlab function to apply (string or function handle)
% arg  - Arguments of the function
%        > Arguments that should be iterated over should be preceded
%          with 'iter'
%        > Arguments that should be iterated in place (i.e. the output
%          replaces the input) should be preceded with 'inplace'
% out  - Output of the function. Each one is a cell array of outputs,
%        unless some arguments were 'inplace'. In this case, the
%        function should return inplace arguments first, in the same
%        order as their input order.
%
% LOCAL MODE
% ----------
% FORMAT [out1, ...] = distribute(func, ('iter'/'inplace'), arg1, ...)
% FORMAT [out1, ...] = distribute(nworkers, func, ('iter'/'inplace'), arg1, ...)
%
% -------------------------------------------------------------------------
% CONFIGURATION
% -------------
% 
% To use this kind of distributed computing, you need a workstation and a
% cluster that have both access to a shared folder. This allows the main
% script and distributed jobs to share data by writing it on disk, rather
% than using copies over SSH.
%
% You also must be able to SSH into the cluster without having to manually
% enter your password, i.e., either the cluster is not password-protected
% or your RSA key must be registered with the cluster. To do this you
% should:
% LINUX
% 1) generate your own RSA key on your workstation
%    >> ssh-keygen -t rsa
% 2) register this key with the cluster
%    >> ssh-copy-id login@cluster
%
% -------------------------------------------------------------------------
% EXAMPLE 1)
% ----------
%
% Let 'a' and 'b' be lists (cell array) of numeric arrays:
% >> a = num2cell(rand(5,5,10),[1 2]);
% >> b = num2cell(rand(5,5,10),[1 2]);
%
% We want to compute all sums a{i} + b{i}. We would call:
% >> [opt,c] = distribute(opt, 'plus', 'iter', a, 'iter', b);
%
% It will perform the equivalent of:
% >> c = cell(size(a));
% >> for i=1:numel(a)
% >>     c{i} = plus(a{i}, b{i});
% >> end
%
% To perform the same operation in place:
% >> [opt,a] = distribute(opt, 'plus', 'inplace', a, 'iter', b);
%
% which is equivalent to:
% >> for i=1:numel(a)
% >>     a{i} = plus(a{i}, b{i});
% >> end
%
% EXAMPLE 2)
% ----------
%
% Let 'dat' be a structure array. We want to apply the function 
% 'processOneData' to each of its elements. Let 'info' be a structure 
% which is useful to all of dat elements. We would call:
% >> [opt,dat] = distribute(opt, 'processOneData', 'inplace', dat, info)
%
% It will perform the equivalent of:
% >> for i=1:numel(dat)
% >>     dat(i) = processOneData(dat(i), info);
% >> end
% _________________________________________________________________________

    if isnumeric(varargin{1}) && isscalar(varargin{1})
        opt = struct;
        opt.client.workers = varargin{1};
        opt.mode = 'parfor';
        opt = distribute_default(opt);
        func = varargin{2};
        varargin  = varargin(3:end);
        no_option = true;
        varargout = cell(1,nargout+1);
    elseif ~isstruct(varargin{1})
        opt  = [];
        func = varargin{1};
        varargin  = varargin(2:end);
        no_option = true;
        varargout = cell(1,nargout+1);
    else
        opt  = varargin{1};
        func = varargin{2};
        varargin  = varargin(3:end);
        no_option = false;
        varargout = cell(1,nargout);
    end
    if isempty(opt)
        opt = struct;
        opt = distribute_default(opt);
    end

    % Parse input
    % -----------
    args  = {};
    flags = {};
    access = {};
    N = 1;
    while ~isempty(varargin)
        if ischar(varargin{1}) && any(strcmpi(varargin{1}, {'iter','inplace'}))
            flags    = [flags {varargin{1}}];
            args     = [args  {varargin{2}}];
            if iscell(varargin{2})
                access = [access  {'{}'}];
            elseif isstruct(varargin{2})
                access = [access {'()'}];
            else
                error(['distribute: an iterable input must either be ' ...
                       'a cell array or a struct array'])
            end
            if N > 1 && numel(varargin{2}) ~= N
                error(['distribute: all iterable inputs should have ' ...
                      'the same number of elements'])
            end
            N = numel(varargin{2});
            varargin = varargin(3:end);
        else
            flags    = [flags {''}];
            args     = [args  {varargin{1}}];
            access   = [access {''}];
            varargin = varargin(2:end);
        end
    end
    
    % Convert function name <-> function handle
    % -----------------------------------------
    if ischar(func)
        funcstr = func; 
        func    = str2func(func);
    else
        funcstr = func2str(func);
    end
    
    
    % Distribute
    % ----------
    if strcmpi(opt.mode, 'qsub') && opt.server.setup && check_server_load(opt)
         
        if ~iscell(opt.job.mem)
            opt.job.mem = {opt.job.mem};
        end
        if ~opt.job.batch && numel(opt.job.mem) < N
            opt.job.mem = padarray(opt.job.mem, [0 N-numel(opt.job.mem)], 'replicate', 'post');
        end
        
        if opt.job.batch
            [varargout{1:numel(varargout)}] = distribute_server_batch(opt, funcstr, args, flags, access, N);
        else
            [varargout{1:numel(varargout)}] = distribute_server_ind(opt, funcstr, args, flags, access, N);
        end
        
        opt = varargout{1};
        
        if opt.job.est_mem || opt.verbose
            % Estimate new memory usage (or just output memory usage)
            % -------------------------
            opt = estimate_mem(opt);
        end
    elseif double(opt.client.workers) == 0 || strcmpi(opt.mode, 'for')
        [varargout{2:numel(varargout)}] = distribute_not(opt, func, args, flags, access, N);
    else
        [varargout{2:numel(varargout)}] = distribute_local(opt, func, args, flags, access, N);
    end

    % Output option structure
    % -----------------------
    if no_option,   varargout    = varargout(2:end);
    else,           varargout{1} = opt;    
    end
end

function ok = check_server_load(~)
    ok = true;
end







