function Z = mrf_post(Z,vs,opt)
% Ad-hoc MRF clean-up of segmentation     
P   = zeros(size(Z),'uint8');   
G   = ones([size(Z,4),1],'single')*opt.clean.mrf.strength;
vx2 = 1./single(vs);
for i=1:opt.clean.mrf.niter        
    spm_mrf(P,Z,G,vx2);        
end       

Z = single(P)/255; 
% figure; imshow3D(squeeze(Z))    
%==========================================================================