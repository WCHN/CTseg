function model = init_template_bg(model,opt)
dm                   = model.template.nii.dat.dim;
K                    = dm(4);
model.template.bg    = cell(1,2);
model.template.bg{1} = zeros(K,1);
model.template.bg{2} = zeros(K,1);

if opt.template.bg_class > 0 && opt.template.bg_class <= K
    mx = max(model.template.nii.dat(:));
    mn = min(model.template.nii.dat(:));
end

for k=1:K
    if opt.template.bg_class > 0 && opt.template.bg_class <= K
        if k == opt.template.bg_class
            model.template.bg{1}(k) = mx;    
            model.template.bg{2}(k) = mx;    
        else
            model.template.bg{1}(k) = mn;    
            model.template.bg{2}(k) = mn;    
        end        
    else
        model.template.bg{1}(k) = mean(mean(model.template.nii.dat(:,:,1,k)));
        model.template.bg{2}(k) = mean(mean(model.template.nii.dat(:,:,end,k)));    
    end
end
% figure;imagesc(model.template.nii.dat(:,:,1,7))
%==========================================================================