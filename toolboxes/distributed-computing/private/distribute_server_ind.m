function varargout = distribute_server_ind(opt, func, args, flags, access, N)
% -------------------------------------------------------------------------
%   Distribute on server (individual)
% -------------------------------------------------------------------------

    % Temporary directory
    % -------------------
    if exist('java.util.UUID', 'class')
        uid = char(java.util.UUID.randomUUID()); % Session UID
    else
        uid = datestr(now, 'yyyymmddTHHMMSS'); % Hopfully no double call in the same second
    end
    client_dir  = fullfile(opt.client.folder, uid);
    server_dir  = fullfile(opt.server.folder, uid);
    mkdir(client_dir);
    
    
    % Replicate memory requirement
    % ----------------------------
    if numel(opt.job.mem) == 1
        [opt.job.mem{1:N}] = deal(opt.job.mem{1});
    end
    
    % Write input arguments
    % ---------------------
    push_arguments(opt, client_dir, args, flags, access, N)
        
    % Write main script
    % -----------------
    main_name = cell(1,N);
    main_sh   = cell(1,N);
    main_end  = cell(1,N);
    for n=1:N
        [main_name{n}, main_sh{n}, main_end{n}] = create_individual_job(opt, client_dir, server_dir, func, n, nargout-1);
    end
    
    % Submit main script
    % ------------------
    opt.job.id = submit_main_job(opt, server_dir, main_sh);
    
    gather_end = '';
    if opt.job.use_dummy
        
        % Write gather script
        % -----------------
        [gather_sh, gather_end] = create_gather_job(opt, client_dir, server_dir);
                           
        % Submit gather job
        % ----------------- 
        dependency = main_name{1};
        for n=2:N
            dependency = [dependency ',' main_name{n}];
        end
        submit_gather_job(opt, server_dir, gather_sh, dependency);
        
    end
    
    if opt.verbose
        uihelper('submitted', N, opt.job.batch, opt.job.id);
    end
    
    % Track jobs
    % ----------
    track_jobs(client_dir, main_end, gather_end, opt.verbose);
    
    % Store opt
    %-----------
    varargout{1} = opt;
    
    % Read output
    % -----------
    [varargout{2:nargout}] = pull_results(opt, client_dir, args, flags, access, N);
    
    % Clean disk
    % ----------
    if opt.clean
        rmdir(client_dir);
    end
    
end