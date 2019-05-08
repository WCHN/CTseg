function Z_a = push_responsibilities(Z,y,dm_a)
K   = size(Z,4);
Z_a = zeros([dm_a K],'single');
for k=1:K
    Z_a(:,:,:,k) = spm_diffeo('push',Z(:,:,:,k),y,dm_a(1:3));
end
%==========================================================================