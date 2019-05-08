function varargout = distribute_local(opt, func, args, flags, access, N)
% -------------------------------------------------------------------------
%   Distribute locally
% -------------------------------------------------------------------------
% It is not very efficient now as all N arguments have to be sent to every
% worker. I could not find a way to be generic AND allow slicing in parfor.

    % Prepare temporary output
    % ------------------------
    Nout = nargout;
    out  = cell(N, Nout);
    
    if opt.verbose
        start_track = uihelper('begin',N);
    end
    
    % Iterate
    % -------
    % /!\ no efficient slicing
    parfor (n=1:N, double(opt.client.workers))
        args1 = cell(size(args));
        for i=1:numel(args)
            switch lower(flags{i})
                case ''
                    args1{i} = args{i};
                case {'iter', 'inplace'}
                    switch lower(access{i})
                        case '()'
                            args1{i} = args{i}(n);
                        case '{}'
                            args1{i} = args{i}{n};
                    end
            end
        end
        out1        = cell(1,Nout);
        [out1{1,:}] = func(args1{:});        
        out(n,:)    = out1;        
    end

    % Write final output
    % ------------------
    [varargout{1:nargout}] = deal({});
    j = 1;
    for i=1:numel(args)
        if strcmpi(flags{i}, 'inplace')
            varargout{j} = args{i};
            if strcmpi(access{i}, '{}')
                for n=1:N
                    varargout{j}{n} = out{n,j};
                end
            elseif strcmpi(access{i}, '()')
                for n=1:N
                    varargout{j}(n) = out{n,j};
                end
            end
            j = j + 1;
        end
    end
    j1 = j;
    for j=j1:nargout
        varargout{j} = out(:,j)';
    end
    
    if opt.verbose
        uihelper('end', N, start_track);
    end
end