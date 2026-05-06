% Convert all .mha files to .nii.gz in the SynthRAD2025 dataset.
% Requires SPM12 on the MATLAB path (for nifti/file_array).

addpath(fileparts(fileparts(mfilename('fullpath'))));  % experiments/
cfg = config();
data_dir = cfg.data_dir;

mha_files = dir(fullfile(data_dir, '**', '*.mha'));
fprintf('Found %d .mha files\n', numel(mha_files));

for i = 1:numel(mha_files)
    mha_path = fullfile(mha_files(i).folder, mha_files(i).name);
    nii_path = strrep(mha_path, '.mha', '.nii.gz');

    if exist(nii_path, 'file')
        fprintf('  SKIP (exists): %s\n', nii_path);
        continue;
    end

    [img, mat, dtype] = read_mha(mha_path);
    write_nii_gz(nii_path, img, mat, dtype);
    fprintf('  %s -> %s\n', mha_files(i).name, strrep(mha_files(i).name, '.mha', '.nii.gz'));
end

fprintf('Done.\n');


function [img, mat, dtype] = read_mha(fname)
% Read a .mha (MetaImage) file and return the image volume,
% a 4x4 voxel-to-world affine matrix, and the output NIfTI datatype.

    fid = fopen(fname, 'r');
    if fid == -1, error('Cannot open %s', fname); end

    % Parse header
    ndims = 3;
    dim   = [1 1 1];
    offset = [0 0 0];
    spacing = [1 1 1];
    R = eye(3);
    elem_type = '';
    header_size = -1;
    is_local = true;

    while true
        line = fgetl(fid);
        if ~ischar(line), break; end

        tokens = strsplit(strtrim(line), '=');
        if numel(tokens) < 2, continue; end
        key = strtrim(tokens{1});
        val = strtrim(strjoin(tokens(2:end), '='));

        switch upper(key)
            case 'NDIMS'
                ndims = str2double(val);
            case 'DIMSIZE'
                dim = sscanf(val, '%d')';
            case 'OFFSET'
                offset = sscanf(val, '%f')';
            case 'ELEMENTSPACING'
                spacing = sscanf(val, '%f')';
            case 'TRANSFORMMATRIX'
                vals = sscanf(val, '%f');
                R = reshape(vals, [ndims ndims])';
            case 'ELEMENTTYPE'
                elem_type = val;
            case 'HEADERSIZE'
                header_size = str2double(val);
            case 'ELEMENTDATAFILE'
                is_local = strcmpi(val, 'LOCAL');
                break;
        end
    end

    % Determine binary data type
    switch upper(elem_type)
        case 'MET_SHORT',  matlab_type = 'int16';   nii_dtype = 'int16';
        case 'MET_USHORT', matlab_type = 'uint16';  nii_dtype = 'int16';
        case 'MET_INT',    matlab_type = 'int32';   nii_dtype = 'int32';
        case 'MET_UINT',   matlab_type = 'uint32';  nii_dtype = 'int32';
        case 'MET_FLOAT',  matlab_type = 'single';  nii_dtype = 'float32';
        case 'MET_DOUBLE', matlab_type = 'double';  nii_dtype = 'float64';
        case 'MET_UCHAR',  matlab_type = 'uint8';   nii_dtype = 'uint8';
        case 'MET_CHAR',   matlab_type = 'int8';    nii_dtype = 'int8';
        otherwise, error('Unsupported ElementType: %s', elem_type);
    end

    if ~is_local
        fclose(fid);
        error('Only LOCAL ElementDataFile is supported');
    end

    % Read binary data from current position (right after header)
    if header_size > 0
        fseek(fid, header_size, 'bof');
    end
    img = fread(fid, prod(dim), ['*' matlab_type]);
    fclose(fid);

    img = reshape(img, dim);
    dtype = nii_dtype;

    % Build voxel-to-world affine (RAS)
    % MHA stores direction cosines in R and spacing separately
    mat = eye(4);
    mat(1:3, 1:3) = R * diag(spacing);
    mat(1:3, 4)   = offset(:);
end


function write_nii_gz(fname, img, mat, dtype)
% Write a NIfTI .nii.gz file using SPM's nifti/file_array.

    % Write uncompressed .nii first, then gzip
    nii_path = regexprep(fname, '\.gz$', '');

    dim = size(img);
    Nii         = nifti;
    Nii.dat     = file_array(nii_path, dim, dtype, 0);
    Nii.mat     = mat;
    Nii.mat0    = mat;
    Nii.descrip = 'Converted from MHA';
    create(Nii);
    Nii.dat(:,:,:) = img;

    gzip(nii_path);
    delete(nii_path);
end
