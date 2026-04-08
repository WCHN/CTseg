function [hd95, assd] = compute_surface_distances(A, B, voxel_size)
% Compute 95th percentile Hausdorff distance and average symmetric surface
% distance between two binary volumes.
%
%   A, B:        3D binary volumes (same size)
%   voxel_size:  [vx vy vz] in mm (default: [1 1 1])
%   hd95:        95th percentile Hausdorff distance (mm)
%   assd:        average symmetric surface distance (mm)
%
% Requires Image Processing Toolbox (bwperim, bwdist).

    if nargin < 3, voxel_size = [1 1 1]; end

    % Handle empty masks
    if ~any(A(:)) || ~any(B(:))
        hd95 = NaN;
        assd = NaN;
        return;
    end

    % Extract surfaces
    surf_A = bwperim(A);
    surf_B = bwperim(B);

    % Distance transform of each volume (in voxels), scaled to mm
    dist_A = bwdist(A);
    dist_B = bwdist(B);

    % Scale distance maps by voxel size for anisotropic voxels.
    % bwdist returns Euclidean distance in voxel units assuming isotropic
    % spacing. We resample the volume to isotropic 1mm spacing before
    % computing distances, then extract surface distances in mm.
    iso_vox = min(voxel_size);  % target isotropic voxel size
    scale = voxel_size / iso_vox;
    if any(abs(scale - 1) > 0.01)
        % Resample to isotropic voxels via nearest-neighbour
        [d1, d2, d3] = size(A);
        new_sz = round([d1 d2 d3] .* scale);
        [X, Y, Z] = ndgrid(linspace(1, d1, new_sz(1)), ...
                            linspace(1, d2, new_sz(2)), ...
                            linspace(1, d3, new_sz(3)));
        A_iso = interp3(double(A), Y, X, Z, 'nearest', 0) > 0.5;
        B_iso = interp3(double(B), Y, X, Z, 'nearest', 0) > 0.5;
        surf_A = bwperim(A_iso);
        surf_B = bwperim(B_iso);
        dist_A = bwdist(A_iso);
        dist_B = bwdist(B_iso);
    end

    % Distances from surface of A to nearest point in B, and vice versa (in mm)
    d_A2B = dist_B(surf_A) * iso_vox;
    d_B2A = dist_A(surf_B) * iso_vox;

    % 95th percentile Hausdorff distance
    hd95 = prctile(double([d_A2B(:); d_B2A(:)]), 95);

    % Average symmetric surface distance
    assd = (mean(double(d_A2B(:))) + mean(double(d_B2A(:)))) / 2;
end
