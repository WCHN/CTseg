function subjects = get_subjects(cfg)
% Return list of subject directory structs.
% Excludes subjects listed in cfg.exclude (incomplete coverage etc.).
% In test mode (cfg.test_mode=true), returns a limited number of subjects.

    d = dir(cfg.data_dir);
    subjects = d([d.isdir] & ~startsWith({d.name}, '.') ...
                            & ~strcmp({d.name}, 'overviews') ...
                            & ~strcmp({d.name}, 'results'));

    % Apply exclusion list
    if isfield(cfg, 'exclude') && ~isempty(cfg.exclude)
        keep = ~ismember({subjects.name}, cfg.exclude);
        n_excluded = sum(~keep);
        subjects = subjects(keep);
        if n_excluded > 0
            fprintf('  Excluded %d subjects (incomplete coverage)\n', n_excluded);
        end
    end

    if cfg.test_mode
        if ~isempty(cfg.test_subjects)
            idx = ismember({subjects.name}, cfg.test_subjects);
            if ~any(idx)
                error('None of the test subjects found in %s', cfg.data_dir);
            end
            subjects = subjects(idx);
        else
            n = min(cfg.test_n_subjects, numel(subjects));
            subjects = subjects(1:n);
        end
        fprintf('  [TEST MODE] Processing %d subject(s): %s\n', ...
            numel(subjects), strjoin({subjects.name}, ', '));
    end
end
