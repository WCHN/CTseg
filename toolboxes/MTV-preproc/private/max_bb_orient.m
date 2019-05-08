function [mat,dm] = max_bb_orient(Nii,vx,padding)
% Calculate orientation matrix and dimensions from maximum bounding-box
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 3, padding = 0; end

mn = [ Inf  Inf  Inf]';
mx = [-Inf -Inf -Inf]';
for i=1:numel(Nii)
    N = numel(Nii{i});
    for n=1:N
        dmn = size(Nii{i}(n).dat);
        
        if numel(dmn) == 2
            dmn(3) = 0; 
        end

        t = uint8(0:7);
        c = diag(dmn+1)*double([bitshift(bitand(t,bitshift(uint8(1),1-1)),1-1)
                              bitshift(bitand(t,bitshift(uint8(1),2-1)),1-2)
                              bitshift(bitand(t,bitshift(uint8(1),3-1)),1-3)]);
        c  = bsxfun(@plus,Nii{i}(n).mat(1:3,1:3)*c,Nii{i}(n).mat(1:3,4));
        mx = max(mx,max(c,[],2));
        mn = min(mn,min(c,[],2));
    end
end

mat = spm_matrix(mn-1)*diag([vx 1])*spm_matrix(-[1 1 1]);
dm  = ceil((mat\[mx'+1 1]')');
dm  = dm(1:3);

if padding > 0
    dm         = dm + 2*padding;
    mat(1:3,4) = mat(1:3,4) - padding;
end

if dmn(3) == 0
    dm(3) = 1;
end
%==========================================================================