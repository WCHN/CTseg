function [Nii,M] = reset_origin(Nii,vx)
if nargin < 2, vx = [];    end

fprintf('Resetting origin...')
N    = numel(Nii{1});
M    = cell(1,N);
M(:) = {eye(4)};
for n=1:N
    f = Nii{1}(n).dat.fname;    
    f = nm_reorient(f,vx);
    
    M{n} = do_reset_origin(f);
    
    Nii{1}(n) = nifti(f);
end

if numel(Nii) > 1
    % Keep labels in alignment
    for n=1:N
        if isempty(Nii{2}(n).dat), continue; end
        
        f = Nii{2}(n).dat.fname;    
        f = nm_reorient(f,vx);
        do_reset_origin(f)

        Nii{2}(n) = nifti(f);
    end    
end
fprintf('done!\n')
%==========================================================================

%==========================================================================
function npth = nm_reorient(pth,vx,prefix,deg)
if nargin < 2, vx     = [];    end
if nargin < 3, prefix = 'ro_'; end
if nargin < 4, deg    = 1;     end

if ~isempty(vx) && length(vx) < 3
    vx=[vx vx vx];
end

% Get information about the image volumes
VV = spm_vol(pth);

for V=VV' % Loop over images

    % The corners of the current volume
    d = V.dim(1:3);
    c = [	1    1    1    1
        1    1    d(3) 1
        1    d(2) 1    1
        1    d(2) d(3) 1
        d(1) 1    1    1
        d(1) 1    d(3) 1
        d(1) d(2) 1    1
        d(1) d(2) d(3) 1]';

    % The corners of the volume in mm space
    tc = V.mat(1:3,1:4)*c;
    if spm_flip_analyze_images, tc(1,:) = -tc(1,:); end

    % Max and min co-ordinates for determining a bounding-box
    mx = round(max(tc,[],2)');
    mn = round(min(tc,[],2)');

    vx0 = sqrt(sum(V.mat(1:3,1:3).^2));
    if isempty(vx)
        vx = vx0;
    else
        vx(vx0 > 1) = vx0(vx0 > 1);
    end    
    
    % Translate so that minimum moves to [1,1,1]
    % This is the key bit for changing voxel sizes,
    % output orientations etc.
    mat = spm_matrix(mn)*diag([vx 1])*spm_matrix(-[1 1 1]);

    % Dimensions in mm
    dim = ceil((mat\[mx 1]')');

    % Output image based on information from the original
    VO               = V;

    % Create a filename for the output image (prefixed by 'r')
    [lpath,name,ext] = fileparts(V.fname);
    VO.fname         = fullfile(lpath,[prefix name ext]);

    % Dimensions of output image
    VO.dim(1:3)      = dim(1:3);

    % Voxel-to-world transform of output image
    if spm_flip_analyze_images, mat = diag([-1 1 1 1])*mat; end
    VO.mat           = mat;

    % Initialise plot of how far reslicing has gone
    %spm_progress_bar('Init',dim(3),'reslicing...','planes completed');

    % Create .hdr and open output .img
    VO = spm_create_vol(VO);

    for i=1:dim(3) % Loop over slices of output image

        % Mapping from slice i of the output image,
        % to voxels of the input image
        M   = inv(spm_matrix([0 0 -i])*inv(VO.mat)*V.mat);

        % Extract this slice according to the mapping
        img = spm_slice_vol(V,M,dim(1:2),deg);

        % Write this slice to output image
        spm_write_plane(VO,img,i);

        % Update the progress bar
        %spm_progress_bar('Set',i);

    end % End loop over output slices

    % Get rid of the progress bar
    %spm_progress_bar('Clear');

end % End loop over images

delete(V.fname);
npth = VO.fname;
%==========================================================================

%==========================================================================
function Mout = do_reset_origin(pth,orig)
if nargin < 2, orig = []; end

V   = spm_vol(pth);
M   = V.mat;
dim = V.dim;
vx  = sqrt(sum(M(1:3,1:3).^2));

if det(M(1:3,1:3))<0
    vx(1) = -vx(1); 
end

if isempty(orig)
    orig = (dim(1:3)+1)/2;
end

off  = -vx.*orig;
M1   = [vx(1) 0      0         off(1)
           0      vx(2) 0      off(2)
           0      0      vx(3) off(3)
           0      0      0      1];

V    = spm_vol(pth);
M0   = V.mat;
Mout = M0/M1;

spm_get_space(pth,M1);   
%==========================================================================