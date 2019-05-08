function opt = estimate_mem(opt)
% -------------------------------------------------------------------------
%   Estimate memory usage
% -------------------------------------------------------------------------

    sd = opt.job.sd;
    
    % ---------------------------------------------------------------------
    % BATCH MODE
    if opt.job.batch                
        jobid = opt.job.id{1};
        omem  = opt.job.mem{1};

        cmd = '';
        for i=1:numel(opt.client.source)
            cmd = [cmd 'source ' opt.client.source{i} ' >/dev/null 2>&1 ; '];
        end
        cmd = [cmd opt.ssh.bin ' ' opt.ssh.opt ' ' opt.server.login '@' opt.server.ip ' "'];
        for i=1:numel(opt.server.source)
            cmd = [cmd 'source ' opt.server.source{i} ' >/dev/null 2>&1 ; '];
        end
        cmd = [cmd opt.sched.acct ' '];
        cmd = [cmd ' -j ' num2str(jobid) ' | grep maxvmem"'];
          
        [status,result] = system(cmd);   

        if status==0
            result = regexp(result, 'maxvmem\W+(?<mem>\d+.?\d*)(?<unit>[KMGT]?)', 'names');
            N = numel(result) - 1;
            
            % Convert all in bytes
            a = zeros(1,N);
            for n=1:N
                a(n) = str2double(result(n).mem);
                switch lower(result(n).unit)
                    case 'k'
                        a(n) = a(n) * 1024;
                    case 'm'
                        a(n) = a(n) * (1024^2);
                    case 'g'
                        a(n) = a(n) * (1024^3);
                    case 't'
                        a(n) = a(n) * (1024^4);
                end
            end
            % Compute maximum value
            [mx, i] = max(a);
            a = (1 + sd)*mx;
            % Use biggest unit
            units   = {result.unit};
            mxunit  = units{i(1)};
            switch lower(mxunit)
                case 'k'
                    a = a ./ 1024;
                case 'm'
                    a = a ./ (1024^2);
                case 'g'
                    a = a ./ (1024^3);
                case 't'
                    a = a ./ (1024^4);
            end
            mem = ceil(a * 10)/10; % Ceil to one decimal place
            
            if opt.job.est_mem
                opt.job.mem{1} = [num2str(mem) mxunit];  
            end
        else
            opt.job.mem{1} = omem; 
        end

        if opt.verbose && ~opt.job.est_mem
            fprintf('Memory usage is %s\n',mem);
        elseif opt.verbose
            fprintf('New memory usage is %s (old memory was %s)\n',...
                    opt.job.mem{1},omem);
        end
        
        
    % ---------------------------------------------------------------------
    % INDIVIDUAL MODE
    else        
        N = numel(opt.job.id);

        cmd = '';
        for i=1:numel(opt.client.source)
            cmd = [cmd 'source ' opt.client.source{i} ' >/dev/null 2>&1 ; '];
        end
        cmd = [cmd opt.ssh.bin ' ' opt.ssh.opt ' ' opt.server.login '@' opt.server.ip ' "'];
        for i=1:numel(opt.server.source)
            cmd = [cmd 'source ' opt.server.source{i} ' >/dev/null 2>&1 ; '];
        end
        for n=1:N
            jobid = opt.job.id{n};
            cmd = [cmd opt.sched.acct ' -j ' jobid ' | grep maxvmem ; '];
        end
        cmd = [cmd '"'];

        [status,result] = system(cmd);   

        if status==0
            result = regexp(result, 'maxvmem\W+(?<mem>\d+.?\d*)(?<unit>[KMGT]?)', 'names');
            for n=1:N
                omem = opt.job.mem{n};
                a   = str2double(result(n).mem);
                a   = (1 + sd)*a;
                mem = ceil(a * 10)/10; % Ceil to one decimal place
                
                if opt.job.est_mem
                    opt.job.mem{n} = [num2str(mem) result(n).unit]; 
                    
                    if opt.verbose
                        fprintf('New memory usage is %s (old memory was %s)\n',...
                            opt.job.mem{n},omem);
                    end
                else
                    fprintf('Memory usage is %s\n',mem);
                end            
            end
        end
    end
end