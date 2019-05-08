function [sh,end_job] = create_gather_job(opt, client_dir, server_dir)

    % ---------------------------------------------------------------------
    % Default names
    name    = 'gather';
    sh      = 'gather.sh';        
    cout    = 'gather_cout.log';  
    cerr    = 'gather_cerr.log'; 
    end_job = 'gather_finished';
    
    % ---------------------------------------------------------------------
    % Script
    bash_script = [ ...
        '#!' opt.sh.bin '\n'...
         '\n'...
         '#$ -S ' opt.sh.bin '\n'...
         '#$ -N ' name '\n'...
         '#$ -o ' fullfile(server_dir, cout) '\n'...
         '#$ -e ' fullfile(server_dir, cerr) '\n'...
         '#$ -j n \n'...
         '#$ -t 1-1 \n'...
         '\n\n'];
              
    % ---------------------------------------------------------------------
    % Finished file
    bash_script = [bash_script ...
        'touch ' fullfile(server_dir, end_job) '\n'];               
    
    % ---------------------------------------------------------------------
    % Write on disk
    pth = fullfile(client_dir,sh);
    fid = fopen(pth,'w');
    fprintf(fid,bash_script);
    fclose(fid);
    fileattrib(pth,'+x','u')
end