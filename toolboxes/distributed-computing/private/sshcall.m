function call = sshcall(opt, cmd, hide)
    if nargin < 3
        hide = false;
    end
    call = [];
    for i=1:numel(opt.client.source)
        call = [call 'source ' opt.client.source{i} ' >/dev/null 2>&1 ; '];
    end
    if ~isempty(opt.server.ip)
        call = [call opt.ssh.bin ' ' opt.ssh.opt ' ' opt.server.login '@' opt.server.ip ' "'];
        for i=1:numel(opt.server.source)
            call = [call 'source ' opt.server.source{i} ' >/dev/null 2>&1; '];
        end
    end
    call = [call cmd];
    if ~isempty(opt.server.ip)
        call = [call '"'];
    end
    if hide
        call = [call ' >/dev/null 2>&1'];
    end
end