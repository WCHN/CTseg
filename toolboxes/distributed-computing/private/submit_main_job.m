function id = submit_main_job(opt, server_dir, sh)

    % ---------------------------------------------------------------------
    % Foormat arguments
    if ~iscell(sh)
        sh = {sh};
    end
    mem = opt.job.mem;
    if ~iscell(mem)
        mem = {mem};
    end
    if numel(mem) == 1
        mem = repmat(mem, 1, numel(sh));
    end

    % ---------------------------------------------------------------------
    % Build command
    cmd = '';
    for i=1:numel(sh)
        cmd = [cmd opt.sched.sub ' '];
        switch lower(opt.sched.type)
            case 'sge'
                cmd = [cmd ' -l vf=' mem{i} ...
                           ' -l h_vmem=' mem{i} ' '];
            otherwise
                error('distribute: scheduler %s not implemented yet', opt.sched.type);
        end
        cmd = [cmd fullfile(server_dir, sh{i}) ' ; '];
    end
    
    % ---------------------------------------------------------------------
    % Call command
    [status, result] = system(sshcall(opt, cmd));
    if status
        fprintf([result '\n']);
        error('distribute: status ~= 0 for main on server!')
    end

    % ---------------------------------------------------------------------
    % Get ID
    s = regexp(result, 'Your job(-array)? (?<id>\d+)', 'names');   
    id = {s.id};
    
end