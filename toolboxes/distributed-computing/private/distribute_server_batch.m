function varargout = distribute_server_batch(opt, func, args, flags, access, N)
% -------------------------------------------------------------------------
%   Distribute on server (batch)
% -------------------------------------------------------------------------
% opt    - Distribute option structre
% func   - Matlab function to apply
% args   - List of arguments
% flags  - Flags for each argument (iter, inplace, ...)
% access - Describe access format for iterable arguments ({}, (), ...)
% N      - Number of jobs

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
    
    % Write input arguments
    % ---------------------
    push_arguments(opt, client_dir, args, flags, access, N)
    
    % Write main script
    % -----------------
    [main_name, main_sh, main_end] = create_batch_job(opt, client_dir, server_dir, func, N, nargout-1);

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
        submit_gather_job(opt, server_dir, gather_sh, main_name);

    end
    
    if opt.verbose
        uihelper('submitted', N, opt.job.batch, opt.job.id);
    end
    
    % Track jobs
    % ----------
    track_jobs(client_dir, main_end, gather_end, opt.verbose);
    
    % Store opt
    % ----------
    varargout{1} = opt;
    
    % Read output
    % -----------
    [varargout{2:nargout}] = pull_results(opt, client_dir, args, flags, access, N);
    
    % Clean disk
    % ----------
    if opt.clean
        rmdir(client_dir, 's');
    end
    
end