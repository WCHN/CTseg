function save_failures(cfg, step_name, failures)
% Save and report failure log for a processing step.
% Appends to cfg.out_dir/failures.mat (accumulates across steps).

    if isempty(failures)
        fprintf('  No failures in %s.\n', step_name);
        return;
    end

    fprintf('\n  ** %s: %d subject(s) failed:\n', step_name, numel(failures));
    for i = 1:numel(failures)
        fprintf('     %s: %s\n', failures{i}.subject, failures{i}.error);
    end
    fprintf('\n');

    % Load existing failures and append
    fail_file = fullfile(cfg.out_dir, 'failures.mat');
    if exist(fail_file, 'file')
        S = load(fail_file);
        all_failures = S.all_failures;
    else
        all_failures = {};
    end
    all_failures = [all_failures, failures];

    save(fail_file, 'all_failures');
    fprintf('  Failure log saved to %s\n', fail_file);
end
