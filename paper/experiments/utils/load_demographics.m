function demo = load_demographics(cfg)
% Load patient demographics (age, sex) from the SynthRAD2025 metadata.
% Returns a struct with fields:
%   demo.ids   - cell array of subject IDs
%   demo.age   - numeric array (NaN where unavailable)
%   demo.sex   - numeric array (1=M, 0=F, NaN=unknown)

    xlsx_file = fullfile(cfg.data_dir, 'overviews', '1_HN_train_parameters.xlsx');
    if ~exist(xlsx_file, 'file')
        error('Metadata file not found: %s', xlsx_file);
    end

    T = readtable(xlsx_file);

    demo.ids = T.ID;
    n = numel(demo.ids);

    % Parse age (format: '073Y' -> 73)
    demo.age = nan(n, 1);
    for i = 1:n
        a = T.PatientAge{i};
        if ~isempty(a) && contains(a, 'Y')
            demo.age(i) = str2double(strrep(a, 'Y', ''));
        end
    end

    % Parse sex (M=1, F=0, other=NaN)
    demo.sex = nan(n, 1);
    for i = 1:n
        s = T.PatientSex{i};
        if strcmp(s, 'M')
            demo.sex(i) = 1;
        elseif strcmp(s, 'F')
            demo.sex(i) = 0;
        end
    end

    fprintf('  Demographics loaded: %d subjects, %d with age, %d with sex\n', ...
        n, sum(~isnan(demo.age)), sum(~isnan(demo.sex)));
end
