function varargout = uihelper(id, varargin)

    switch lower(id)
        case 'begin'
            [varargout{1:nargout}] = batch_begin(varargin{:});
        case 'incr'
            [varargout{1:nargout}] = batch_incr(varargin{:});
        case 'end'
            [varargout{1:nargout}] = batch_end(varargin{:});
        case 'submitted'
            [varargout{1:nargout}] = submitted(varargin{:});
        case 'pulled'
            [varargout{1:nargout}] = pulled(varargin{:});
        case 'sec2ydhms'
            [varargout{1:nargout}] = sec2ydhms(varargin{:});
    end

end

function start = batch_begin(total)
    fprintf(1, '%s | ', datestr(now,'mmmm dd, yyyy HH:MM:SS'))
    n    = floor(log10(total)) + 1;
    s    = sprintf(['%' num2str(n) 'd of %' num2str(n) 'd jobs processed'], 0, total);
    fprintf(1, ['%-' num2str(2*n + 20) 's'], s);
    start = tic;
end

function batch_incr(cur,total)
% cur    - Current number of processed subjects
% total  - Total number of subjects

    n = floor(log10(total)) + 1;
    fprintf(1, repmat('\b',1,2*n + 20));
    s = sprintf(['%' num2str(n) 'd of %' num2str(n) 'd jobs processed'], cur, total);
    fprintf(1, ['%-' num2str(2*n + 20) 's'], s);
end


function batch_end(total, start)
    n = floor(log10(total)) + 1;
    dur = sec2ydhms(toc(start), true);
    fprintf(1, repmat('\b',1,2*n + 20));
    s = sprintf(['%' num2str(n) 'd of %' num2str(n) 'd jobs processed in %s'], total, total, dur);
    fprintf(['%-' num2str(2*n + 20) 's\n'], s);
end

function submitted(N, batch, jobid)            
    date = datestr(now,'mmmm dd, yyyy HH:MM:SS');
    if batch
        fprintf('%s | Batch job submitted to cluster (N = %i, id = %s)\n', date, N, jobid{1});
    else
        fprintf('%s | Individual jobs submitted to cluster (N = %i, id = %s - %s)\n', date, N, jobid{1}, jobid{end});
    end
end

function pulled(dur)
    dur = sec2ydhms(dur, true);
    date = datestr(now,'mmmm dd, yyyy HH:MM:SS');
    fprintf([sprintf('%s | Data pulled in ', date) dur '\n']);
end

function str_end_2 = sec2ydhms(time, compact)
% FORMAT  str = sec2ydhms(duration, compact)
% duration - duration in seconds
% compact  - Use compact format [false]
%
% Convert a duration in seconds (obtained from tic/toc) to a character
% representation in terms of years, dat, hour, minute, seconds.
%
% Default representation: 3 hours, 15 minutes, 1 second
% Compact representation: 3h 15m 1s

    if nargin < 2
        compact = false;
    end
    
    dur = duration(0,0,time);
    elapsed = floor(years(dur));
    dur = dur - years(elapsed(end));
    elapsed = [elapsed floor(days(dur))];
    dur = dur - days(elapsed(end));
    elapsed = [elapsed floor(hours(dur))];
    dur = dur - hours(elapsed(end));
    elapsed = [elapsed floor(minutes(dur))];
    dur = dur - minutes(elapsed(end));
    elapsed = [elapsed floor(seconds(dur))];
    if compact
        units = {'y' 'd' 'h' 'm' 's'};
        space = '';
    else
        units = {'year' 'day' 'hour' 'minute' 'second'};
        space = ' ';
    end
    str_end_2 = '';
    for i=1:numel(elapsed)
        if elapsed(i) > 0
            str_end_2 = [str_end_2 sprintf('%d%s%s', elapsed(i), space, units{i})];
            if ~compact && elapsed(i) > 1
                str_end_2 = [str_end_2 's'];
            end
            if ~compact && sum(elapsed(i+1:end)) > 0
                str_end_2 = [str_end_2 ', '];
            end
        end
    end
    if sum(elapsed) == 0
        str_end_2 = [str_end_2 sprintf('< 0%s%s', space, units{5})];
    end

end