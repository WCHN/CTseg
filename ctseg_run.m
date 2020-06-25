function varargout = ctseg_run(job, action)

if nargin == 1, action = 'run'; end

switch lower(action)
    case 'run'
        varargout{1} = run_job(job);
    case 'vout'
        varargout{1} = vout_job(job);
    otherwise
        error('Unknown argument ("%s").', action);
end
%==========================================================================

%==========================================================================
% Run
%==========================================================================
function vout = run_job(job)

vout = vout_job(job);

if isempty(job.odir{1})  
    odir = '';   
else
    odir = job.odir{1}; 
end
if isempty(job.tc) 
    tc = true(3, 3); 
else
    tc = job.tc;    
end
if isempty(job.def) 
    def = true; 
else
    def = job.def; 
end
if isempty(job.correct_header) 
    correct_header = false; 
else
    correct_header = job.correct_header; 
end

N       = size(job.data,1);
results = cell(1,N);
for n=1:N
    in         = deblank(job.data{n});    
    results{n} = spm_segment_ct(in, odir, tc, def, correct_header);
end
%==========================================================================

%==========================================================================
% Vout
%==========================================================================
function vout = vout_job(job)

if isempty(job.correct_header) 
    correct_header = false; 
else
    correct_header = job.correct_header; 
end

n     = numel(job.data);
parts = cell(n,4);
for j=1:n
    [parts{j,:}] = spm_fileparts(job.data{j});
end

tiss = struct('c',{},'wc',{},'mwc',{});
for i=1:3
    if job.tc(1)
        tiss(i).c = cell(n,1);
        for j=1:n
            if correct_header
                tiss(i).c{j} = fullfile(parts{j,1},['c',num2str(i),'r',parts{j,2},'.nii']);
            else
                tiss(i).c{j} = fullfile(parts{j,1},['c',num2str(i),parts{j,2},'.nii']);
            end
        end
    end
    if job.tc(2)
        tiss(i).wc = cell(n,1);
        for j=1:n
            if correct_header
                tiss(i).wc{j} = fullfile(parts{j,1},['wc',num2str(i),'r',parts{j,2},'.nii']);
            else
                tiss(i).wc{j} = fullfile(parts{j,1},['wc',num2str(i),parts{j,2},'.nii']);
            end
        end
    end
    if job.tc(3)
        tiss(i).mwc = cell(n,1);
        for j=1:n
            if correct_header
                tiss(i).mwc{j} = fullfile(parts{j,1},['mwc',num2str(i),'r',parts{j,2},'.nii']);
            else
                tiss(i).mwc{j} = fullfile(parts{j,1},['mwc',num2str(i),parts{j,2},'.nii']);
            end
        end
    end
end

if job.def
    fordef = cell(n,1);
    for j=1:n
        fordef{j} = fullfile(parts{j,1},['y_',parts{j,2},'.nii']);
    end
else
    fordef = {};
end

vout  = struct('tiss',tiss,'fordef',{fordef});