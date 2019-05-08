function id = submit_gather_job(opt, server_dir, sh, dependency)

    % ---------------------------------------------------------------------
    % Build command
    cmd = [opt.sched.sub ' '];
    cmd = [cmd ...
            ' -l vf=0.1G -l h_vmem=0.1G' ...
            ' -hold_jid ' dependency ...
            ' -cwd '      fullfile(server_dir, sh)];

    % ---------------------------------------------------------------------
    % Call command
    [status,result] = system(sshcall(opt, cmd));
    if status
        fprintf([result '\n'])
        error('status~=0 for gathering job on server!') 
    end

    % ---------------------------------------------------------------------
    % Get ID
    s  = regexp(result, '^\D*(?<id>\d+)', 'names'); % ID of array job on server    
    id = s.id;
end