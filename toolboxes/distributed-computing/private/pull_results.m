function varargout = pull_results(opt, client_dir, args, flags, access, N)

    % Reverse translation
    opt.server_to_client = true;
    
    % initialise output structure
    [varargout{1:nargout}] = deal({});
    j = 1;
    for i=1:numel(args)
        if strcmpi(flags{i}, 'inplace')
            varargout{j} = args{i};
            j = j + 1;
        end
    end
    
    % fill output structure
    for n=1:N
        fargout = fullfile(client_dir, sprintf('out_%d.mat', n));
        % read argout
        try
            load(fargout, 'argout');
            argout = distribute_translate(opt, argout); 
        catch ME
            warning('Error reading file %d (%s)\n', n, fargout)
            for i=1:numel(ME.stack)
                disp([ME.stack(i).name ', line ' num2str(ME.stack(i).line)]);
            end
            disp(ME.message)  
            continue
        end
        % fill inplace
        j = 1;
        for i=1:numel(args)
            if strcmpi(flags{i}, 'inplace')
                if strcmpi(access{i}, '{}')
                    varargout{j}{n} = argout{j};
                elseif strcmpi(access{i}, '()')
                    varargout{j}(n) = argout{j};
                end
                j = j + 1;
            end
        end
        % fill remaining
        j1 = j;
        for j=j1:nargout
            varargout{j}{n} = argout{j};
        end
        clear argout
    end
    
end