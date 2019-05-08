function [name, sh, finished] = create_individual_job(opt, client_dir, server_dir, func, n, Nout)
% [name, all_finished] = create_batch_job(opt, dir, argfile, N)
% opt     - Distribute option structure
% dir     - Directory for temporary data
% argfile - Matlab file containing all input argument path

    % ---------------------------------------------------------------------
    % Filenames
    name       = ['main_' num2str(n)];
    sh         = [name '.sh'];                 % main bash script
    cout       = [name '_cout.log'];           % main output file
    cerr       = [name '_cerr.log'];           % main error file
    argin      = ['in_'  num2str(n) '.mat'];
    argout     = ['out_' num2str(n) '.mat'];
    finished   = [name '_finished'];

    % ---------------------------------------------------------------------
    % SGE header
    script = [             ...
        '#!' opt.sh.bin '\n'     ...
        '\n'                     ...
        '#$ -S ' opt.sh.bin '\n' ... % Shell path
        '#$ -N ' name '\n'       ... % Job name
        '#$ -o ' fullfile(server_dir,cout) '\n'  ... % Path to output file
        '#$ -e ' fullfile(server_dir,cerr) '\n'  ... % Path to error file
        '#$ -j n \n'                ... % Do not join out/err files
        '\n'];
    
    % ---------------------------------------------------------------------
    % Matlab command
    script = [script 'matlab_cmd="'];
    if ~isempty(opt.matlab.priv.add)
        script = [script 'addpath(' opt.matlab.priv.add ');'];
    end
    if ~isempty(opt.matlab.priv.addsub)
        script = [script 'addpath(' opt.matlab.priv.addsub ');'];
    end
    script = [script ...
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
    script = [script ...
        'touch ' fullfile(server_dir, finished) '\n'];
    
    % ---------------------------------------------------------------------
    % Write on disk
    path = fullfile(client_dir, sh);
    fid = fopen(path, 'w');
    fprintf(fid, script);
    fclose(fid);
    fileattrib(path, '+x', 'u');

end