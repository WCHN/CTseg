function push_arguments(opt, dir, args, flags, access, N)
% FORMAT write_input_arg(opt, dir, args, flags, access)
%
% opt    - Distribute option structure
% dir    - Directory were to write data (client side)
% args   - Arguments
% flags  - Flags about arguments (iter, inplace)
% access - Arguments access type ({}, ())
% N      - Number of subjects
%
% Write arguments for each subject in 'in_%d.mat'

    inprefix = 'in_';
    for n=1:N
        
        argin = cell(size(args));
        for i=1:numel(args)
            switch lower(flags{i})
                case ''
                    argin{i} = args{i};
                case {'iter', 'inplace'}
                    switch lower(access{i})
                        case '()'
                            argin{i} = args{i}(n);
                        case '{}'
                            argin{i} = args{i}{n};
                    end
            end
            
        end
        argin = distribute_translate(opt, argin);
        save(fullfile(dir, [inprefix  num2str(n) '.mat']), 'argin', '-mat'); 
        clear argin
    end
end