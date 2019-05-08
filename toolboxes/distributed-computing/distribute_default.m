function opt = distribute_default(opt)
% _________________________________________________________________________
%
%           Default options for distributing Matlab jobs
%
% -------------------------------------------------------------------------
%
% FORMAT function opt = distribute_default(opt)
%
% opt is a structure with optional fields:
%
% GENERAL
% -------
% mode    - Parallelisation mode: 'qsub'/'parfor'/['for']
% verbose - Speak during processing [false]
%
% CLUSTER
% -------
% server.ip     - IP adress (or alias name) of the cluster ['' = cluster == client]
% server.login  - Login with which to connect ['']
% server.source - Files to source on server side [auto]
%                 > Try to find bashrc and/or bash_profile
% server.folder - Shared folder for writing data, scripts, ...
%
% LOCAL
% -----
% client.source  - Files to source on server side [auto]
% client.workers - Number of local workers [auto]
% client.folder  - Shared folder for writing data, scripts, ...
%
% SUBMIT JOBS
% -----------
% ssh.type      - SSH software to use 'ssh'/'putty'/[try to detect]
% ssh.bin       - Path to the ssh binary [try to detect]
% ssh.opt       - SSH options ['-x']
% sched.sub     - Path to the submit binary [try to detect]
% sched.stat    - Path to the stat binary [try to detect]
% sched.acct    - Path to the acct binary [try to detect]
% sched.type    - Type of scheduler 'sge'/'pbs'/[try to detect]
% job.batch     - Submit jobs as a batch (force same mem for all) [true]
% job.mem       - (Initial) Max memory usage by a single job ['2G']
% job.est_mem   - Estimate max memory usage from previous runs [true]
% job.sd        - Amount of extra memory to add to estimated max memory [0.1]
% job.use_dummy - Uses a dummy job to decide when job have finished [false]
%
% MATLAB
% ------
% matlab.bin    - Path to matlab binary [try to detect]
% matlab.add    - Paths to add to Matlab path [{}]
% matlab.addsub - Paths to add to Matlab path, with subdirectories [{}]
% matlab.opt    - Commandline options to pass to matlab
%                 ['-nojvm' '-nodesktop' '-nosplash' '-singleCompThread']
% spm.path      - Path to SPM []
% spm.toolboxes - List of SPM toolboxes to add to Matlab path [{}]
%
% DATA
% ----
% translate - Cell array of size 2xN with translation between client and
%             server paths [{client.folder server.folder}].
%             Example:
%                  {'/home/me/'     '/mnt/users/me/' ;
%                   '/shared/data/' '/mnt/shared/data'}
% restrict   - Restrict translation to a class: 'char'/'file_array'/['']
% clean      - Clean tmp data when finished [true]
% clean_init - Initially clean tmp data [false]
% _________________________________________________________________________

    if nargin < 1
        opt = struct;
    end

    if ~isfield(opt, 'mode')
        opt.mode = 'for';
    end
    if ~isfield(opt, 'verbose')
        opt.verbose = false;
    end
    
    % CLUSTER
    % -------
    opt.server.setup = true;
    if ~isfield(opt, 'server')
        opt.server = struct;
    end
    if ~isfield(opt.server, 'ip')
        opt.server.ip = '';
    end
    if ~isfield(opt.server, 'login')
        opt.server.login = '';
    end
    if ~isfield(opt.server, 'folder')
        opt.server.folder = '~/.distribute';
    end
    if ~isfield(opt, 'client')
        opt.client = struct;
    end
    if ~isfield(opt.client, 'source')
        opt.client.source = {};
        if isunix
            if exist('~/.bash_profile', 'file')
                opt.client.source = [opt.client.source {'~/.bash_profile'}];
            end
            if exist('~/.bashrc', 'file')
                opt.client.source = [opt.client.source {'~/.bashrc'}];
            end
            if exist('/etc/profile', 'file')
                opt.client.source = [opt.client.source {'/etc/profile'}];
            end
        end
    end
    if ~iscell(opt.client.source)
        opt.client.source = {opt.client.source};
    end
    if ~isfield(opt.client, 'workers')
        if strcmpi(opt.mode, 'parfor')
            myCluster          = parcluster;
            opt.client.workers = myCluster.NumWorkers;
            clear myCluster
        else
            opt.client.workers = 0;
        end
    end
    if ~isfield(opt.client, 'folder')
        opt.client.folder = opt.server.folder;
    end

    
    % IF LOCAL MODE: NO NEED TO CONTINUE
    % ----------------------------------
    if any(strcmpi(opt.mode, {'for','parfor'}))
        return
    end

    % SUBMIT JOBS
    % -----------
    if ~isfield(opt, 'ssh')
        opt.ssh = struct;
    end
    if ~isfield(opt.ssh, 'type') && ~isfield(opt.ssh, 'bin')
        [opt.ssh.bin, opt.ssh.type] = auto_detect('ssh', opt);
    elseif isfield(opt.ssh, 'type')
        opt.ssh.bin = auto_detect('ssh', opt.ssh.type, opt);
    elseif isfield(opt.ssh, 'bin')
        [~, name, ~] = fileparts(opt.ssh.bin);
        if any(strcmpi(name, {'putty', 'ssh'}))
            opt.ssh.type = lower(name);
        else
            warning('Cannot detect ssh type\n')
        end
    end
    if isempty(opt.ssh.bin) && ~isempty(opt.server.ip)
        warning('Could not find an ssh binary')
        opt.server.setup = false;
    end
    if ~isfield(opt.ssh, 'opt')
        if strcmpi(opt.ssh.type, 'ssh')
            opt.ssh.opt = '-x';
        else
            opt.ssh.opt = '';
        end
    end
    if ~isfield(opt.server, 'source')
        % Need ssh for that
        opt.server.source = auto_detect('source', opt);
    end
    if ~iscell(opt.server.source)
        opt.server.source = {opt.server.source};
    end
    if ~isfield(opt, 'sched')
        opt.sched = struct;
    end
    if ~isfield(opt.sched, 'sub') ...
            || ~isfield(opt.sched, 'stat') ...
            || ~isfield(opt.sched, 'acct')
        [opt.sched.sub, opt.sched.stat, opt.sched.acct] = auto_detect('sched', opt);
        if isempty(opt.sched.sub)
            opt.server.setup = false;
        end
    end
    if ~isfield(opt.sched, 'type')
        opt.sched.type = auto_detect('sched', 'type', opt.sched.sub, opt);
    end
    if ~isfield(opt, 'job')
        opt.job = struct;
    end
    if ~isfield(opt.job, 'mem')
        opt.job.mem = {'2G'};
    end
    if ~iscell(opt.job.mem)
        opt.job.mem = {opt.job.mem};
    end
    if ~isfield(opt.job, 'batch')
        opt.job.batch = true;
    end
    if ~isfield(opt.job, 'est_mem')
        opt.job.est_mem = true;
    end
    if ~isfield(opt.job, 'sd')
        opt.job.sd = 0.1;
    end
    if ~isfield(opt.job, 'use_dummy')
        opt.job.use_dummy = false;
    end
    if ~isfield(opt, 'optim')
        opt.optim = struct;
    end
    if ~isfield(opt.optim, 'optim')
        opt.optim.optim = true;
    end
    if ~isfield(opt.optim, 'busy')
        opt.optim.busy = 0.9;
    end
    if ~isfield(opt, 'sh')
        opt.sh = struct;
    end
    if ~isfield(opt.sh, 'bin')
        opt.sh.bin = '/bin/sh'; % I should actually try to detect it
    end

    % MATLAB
    % ------
    if ~isfield(opt, 'matlab')
        opt.matlab = struct;
    end
    if ~isfield(opt.matlab, 'bin')
        opt.matlab.bin = auto_detect('matlab', opt);
        if isempty(opt.matlab.bin)
            opt.server.setup = false;
        end
    end
    if ~isfield(opt.matlab, 'add')
        opt.matlab.add = {};
    end
    if ~iscell(opt.matlab.add)
        opt.matlab.add = {opt.matlab.add};
    end
    if ~isfield(opt.matlab, 'addsub')
        opt.matlab.addsub = {};
    end
    if ~iscell(opt.matlab.addsub)
        opt.matlab.addsub = {opt.matlab.addsub};
    end
    if ~isfield(opt.matlab, 'opt')
        opt.matlab.opt = '-nojvm -nodesktop -nosplash -singleCompThread';
    end
    if ~isfield(opt, 'spm')
        opt.spm = struct;
    end
    if ~isfield(opt.spm, 'path')
        opt.spm.path = auto_detect('spm', opt);
    end
    if ~isfield(opt.spm, 'toolboxes')
        opt.spm.toolboxes = {};
    end
    if ~iscell(opt.spm.toolboxes)
        opt.spm.toolboxes = {opt.spm.toolboxes};
    end


    % DATA
    % ----
    if ~isfield(opt, 'translate')
        opt.translate = {};
    end
    if ~isfield(opt, 'restrict')
        opt.restrict = '';
    end
    if ~isfield(opt, 'clean')
        opt.clean = true;
    end
    if ~isfield(opt, 'clean_init')
        opt.clean_init = false;
    end
    
    % BUILD ADDPATH STRING
    % --------------------
    opt.matlab.priv.add = '';
    for i=1:numel(opt.matlab.add)
        opt.matlab.priv.add = [opt.matlab.priv.add '''' opt.matlab.add{i} ''','];
    end
    if ~isempty(opt.spm.path)
        opt.matlab.priv.add = [opt.matlab.priv.add '''' opt.spm.path ''','];
        for i=1:numel(opt.spm.toolboxes)
            opt.matlab.priv.add = [opt.matlab.priv.add ...
            ''''  fullfile(fullfile(opt.spm.path, 'toolbox'), opt.spm.toolboxes{i}) ''','];
        end
    end
    if ~isempty(opt.matlab.priv.add)
        opt.matlab.priv.add = opt.matlab.priv.add(1:end-1);
    end
    opt.matlab.priv.addsub = '';
    for i=1:numel(opt.matlab.addsub)
        opt.matlab.priv.addsub = [opt.matlab.priv.addsub 'genpath(''' opt.matlab.addsub{i} '''),'];
    end
    if ~isempty(opt.matlab.priv.addsub)
        opt.matlab.priv.addsub = opt.matlab.priv.addsub(1:end-1);
    end

    if ~isempty(opt.server.ip)
        if opt.clean_init && exist(opt.client.folder,'dir')
            rmdir(opt.client.folder,'s');
            mkdir(opt.client.folder);
        elseif ~exist(opt.client.folder,'dir')
            mkdir(opt.client.folder);
        end
    end
end

% =========================================================================
%     AUTO DETECT
% =========================================================================

% -------------------------------------------------------------------------
%   general
% -------------------------------------------------------------------------

function varargout = auto_detect(id, varargin)
    switch lower(id)
        case 'ssh'
            [varargout{1:nargout}] = auto_detect_ssh(varargin{:});
        case 'source'
            [varargout{1:nargout}] = auto_detect_source(varargin{:});
        case 'sched'
            [varargout{1:nargout}] = auto_detect_sched(varargin{:});
        case 'matlab'
            [varargout{1:nargout}] = auto_detect_matlab(varargin{:});
        case 'spm'
            [varargout{1:nargout}] = auto_detect_spm(varargin{:});
        otherwise
            warning('auto_detect: unknown action %s\n', id);
    end
end

function ok = sshexist(opt, file)
    st = system(sshcall(opt, ['find ''' file ''''], true));
    ok = (st == 0);
end

function path = sshpath(opt)
    [~, path] = system(sshcall(opt, 'echo \$PATH"'));
    path = split(deblank(path), newline);
    path = split(path{end}, ':');
end

function ok = sshcommandst(opt, cmd)
    st = system(sshcall(opt, cmd, true));
    ok = (st == 0);
end

function path = sshwhich(opt, bin)
    [st, path] = system(sshcall(opt, ['which ' bin]));
    if st
        path = '';
    else
        path = deblank(path);
    end
end

% -------------------------------------------------------------------------
%   ssh
% -------------------------------------------------------------------------

function [bin, type] = auto_detect_ssh(type, opt)
    if nargin == 1
        opt = type;
        type = '';
    end
    if isempty(opt.server.ip)
        % No cluster so no need for ssh
        bin = '';
        return
    end
    bin = '';
    if nargin == 2
        % Known type
        if ispc
            if exist('C:\', 'dir')
                bin = 'C:\';
            else
                bin = '';
                return
            end
            if exist(fullfile(bin, 'Program Files'), 'dir')
                bin = fullfile(bin, 'Program Files');
            elseif exist(fullfile(bin, 'Program Files (x86)'), 'dir')
                bin = fullfile(bin, 'Program Files (x86)');
            else
                bin = '';
                return
            end
            if exist(fullfile(bin, type), 'dir')
                bin = fullfile(bin, type);
            else
                bin = '';
                return
            end
            if exist(fullfile(bin, [type '.exe']), 'file')
                bin = fullfile(bin, [type '.exe']);
            else
                bin = '';
                return
            end

        elseif isunix
            call = '';
            for j=1:numel(opt.client.source)
                call = [call 'source ' opt.client.source{j} ' >/dev/null 2>&1 ; '];
            end
            call = [call 'echo $PATH'];
            [~, path] = system(call);
            path = split(deblank(path), newline);
            path = split(path{end}, ':');
            for i=1:numel(path)
                if exist(fullfile(path{i}, type), 'file')
                    bin = fullfile(path{i}, type);
                    return
                end
            end
        end
    else
        % Unknown type
        types = {'ssh', 'putty'};
        for i=1:numel(types)
            [bin, type] = auto_detect_ssh(types{i}, opt);
            if ~isempty(bin)
                return
            end
        end
        bin  = '';
        type = '';
    end
end

% -------------------------------------------------------------------------
%   source
% -------------------------------------------------------------------------

function path = auto_detect_source(~)
    path = {};
    if isunix
        if exist('~/.bash_profile', 'file')
            path = [path {'~/.bash_profile'}];
        end
        if exist('~/.bashrc', 'file')
            path = [path {'~/.bashrc'}];
        end
        if exist('/etc/profile', 'file')
            path = [path {'/etc/profile'}];
        end
    end

end

% -------------------------------------------------------------------------
%   sched
% -------------------------------------------------------------------------

function varargout = auto_detect_sched(varargin)

    if nargin == 1
        % find sub/stat/acct
        opt = varargin{1};
        sub  = sshwhich(opt, 'qsub');
        stat = sshwhich(opt, 'qstat');
        acct = sshwhich(opt, 'qacct');
        if isempty(sub) || isempty(stat) || isempty(acct)
            path = sshpath(opt);
            for i=1:numel(path)
                if isempty(sub) && sshexist(opt, [path{i} '/qsub'])
                    sub = [path{i} '/qsub'];
                end
                if isempty(stat) && sshexist(opt, [path{i} '/qstat'])
                    stat = [path{i} '/qstat'];
                end
                if isempty(acct) && sshexist(opt, [path{i} '/qacct'])
                    acct = [path{i} '/qacct'];
                end
                if ~isempty(sub) && ~isempty(stat) && ~isempty(acct)
                    break
                end
            end
        end
        if isempty(sub) || isempty(stat) || isempty(acct)
            warning('Could not detect qsub/qstat/qacct')
        end
        varargout{1} = sub;
        varargout{2} = stat;
        varargout{3} = acct;
    elseif nargin == 3
        % find type
        sub = varargin{2};
        opt = varargin{3};
        type = '';
        if contains(sub, 'gridengine')
            type = 'sge';
        elseif sshcommandst(opt, 'man qsub | grep -i ''pbs''')
            type = 'pbs';
        elseif sshcommandst(opt, 'man qsub | grep -i ''sge''')
            type = 'sge';
        else
            warning(['Could not detect sheduler type for sure. ' ...
                     'Trying SGE.']);
        end
        varargout{1} = type;
    end

end

% -------------------------------------------------------------------------
%   matlab
% -------------------------------------------------------------------------

function bin = auto_detect_matlab(opt)
    bin = '';
    path = sshpath(opt);
    for i=1:numel(path)
        if sshexist(opt, [path{i} '/matlab'])
            bin = [path{i} '/matlab'];
            return
        end
    end
    warning('Cannot detect matlab')
end

% -------------------------------------------------------------------------
%   spm
% -------------------------------------------------------------------------

function path = auto_detect_spm(opt)
    path = '';
end

% -------------------------------------------------------------------------
%   Backward compatibility
% -------------------------------------------------------------------------
function out = newline
    if exist('newline','builtin')
        out = builtin('newline');
        return
    end
    out = char(10);    
end
