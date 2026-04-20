function pth = resolve_atlas(mu)
% Resolve atlas shorthand name to a file path.
% If mu is already a valid file path, returns it unchanged.
% If mu is a shorthand (e.g., 'spm15'), resolves to models/ directory.
% Downloads the atlas if not already present.
%
% FORMAT pth = resolve_atlas(mu)

if exist(mu, 'file')
    pth = mu;
    return
end

dir_models = fullfile(fileparts(which('spm_CTseg')), 'models');
if strcmp(mu, 'ctseg')
    pth = fullfile(dir_models, 'mu_CTseg.nii');
else
    pth = fullfile(dir_models, ['mu_CTseg_' mu '.nii']);
end

if ~exist(pth, 'file')
    % Trigger download by calling spm_CTseg's download mechanism
    fprintf('Atlas ''%s'' not found locally, triggering download...\n', mu);
    % Quick call that downloads but doesn't run segmentation
    try
        spm_CTseg(pth, '', false(6,3), false, false, false, NaN, [], [], mu);
    catch
        % Expected to fail (invalid input), but atlas should be downloaded
    end
    if ~exist(pth, 'file')
        error('Could not resolve atlas ''%s''. File not found: %s', mu, pth);
    end
end
