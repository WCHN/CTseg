function [name, sh, finished] = create_batch_job(opt, client_dir, server_dir, func, N, Nout)
% [name, all_finished] = create_batch_job(opt, dir, argfile, N)
% opt     - Distribute option structure
% dir     - Directory for temporary data
% argfile - Matlab file containing all input argument path

    % ---------------------------------------------------------------------
    % Filenames
    name       = 'main';
    sh         = 'main.sh';                 % main bash script
    cout       = 'main_cout.log';           % main output file
    cerr       = 'main_cerr.log';           % main error file
    argin      = 'in_$SGE_TASK_ID.mat';
    argout     = 'out_$SGE_TASK_ID.mat';
    end_prefix = 'main_finished_';
    end_job    = [end_prefix '$SGE_TASK_ID'];
    finished = cell(1,N);
    for n=1:N
        finished{n} = [end_prefix num2str(n)];
    end
    

    % ---------------------------------------------------------------------
    % SGE header
    batch_script = [             ...
        '#!' opt.sh.bin '\n'     ...
        '\n'                     ...
        '#$ -S ' opt.sh.bin '\n' ... % Shell path
        '#$ -N ' name '\n'       ... % Job name
        '#$ -o ' fullfile(server_dir,cout) '\n'  ... % Path to output file
        '#$ -e ' fullfile(server_dir,cerr) '\n'  ... % Path to error file
        '#$ -j n \n'                ... % Do not join out/err files
        '#$ -t 1-' num2str(N) ' \n' ... % Number of subjobs
        '\n'];
    
    % ---------------------------------------------------------------------
    % Matlab command
    batch_script = [batch_script 'matlab_cmd="'];
    if ~isempty(opt.matlab.priv.add)
        batch_script = [batch_script 'addpath(' opt.matlab.priv.add ');'];
    end
    if ~isempty(opt.matlab.priv.addsub)
        batch_script = [batch_script 'addpath(' opt.matlab.priv.addsub ');'];
    end
    batch_script = [batch_script ...
            'load(fullfile(''' server_dir ''',''' argin '''),''argin'');' ...
            'argout=cell(1,' num2str(Nout) ');' ...
            'func=str2func(''' func ''');' ...
            '[argout{1:' num2str(Nout) '}]=func(argin{:});' ...
            'save(fullfile(''' server_dir ''',''' argout '''),''argout'',''-mat'');' ...
            'quit;' ...
        '"\n' ...
        opt.matlab.bin ' ' opt.matlab.opt ' -r $matlab_cmd \n' ... 
    ];
    
    % ---------------------------------------------------------------------
    % Finished file
    batch_script = [batch_script ...
        'touch ' fullfile(server_dir, end_job) '\n'];
    
    % ---------------------------------------------------------------------
    % Write on disk
    path = fullfile(client_dir, sh);
    fid = fopen(path, 'w');
    fprintf(fid, batch_script);
    fclose(fid);
    fileattrib(path, '+x', 'u');

end