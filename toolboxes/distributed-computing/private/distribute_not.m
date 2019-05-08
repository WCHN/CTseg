function varargout = distribute_not(opt, func, args, flags, access, N)
% -------------------------------------------------------------------------
%   Do not distribute
% -------------------------------------------------------------------------

    % Prepare temporary output
    % ------------------------
    out = cell(N, nargout);
    
    if opt.verbose
        start_track = uihelper('begin',N);
    end
    
    % Iterate
    % -------
    for n=1:N
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
        [out{n,:}] = func(args1{:});
        if opt.verbose
            uihelper('incr', n, N);
        end
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