function [g,H] = mnom_objfun_slice(Z,Template,y1,z,log_prop)
dm   = size(Z);
emu  = cell(dm(4),1);
mu1  = cell(dm(4),1);
semu = 0;
dmu  = cell(dm(4),3);
for k=1:dm(4)
    [mu1{k},dmu{k,1},dmu{k,2},dmu{k,3}] = spm_diffeo('bsplins',Template(:,:,:,k),y1(:,:,z,:),[1 1 1  0 0 0]);
        
    emu{k} = exp(mu1{k} + log_prop(k));
    semu   = semu       + emu{k};
end

msk = sum(Z,4) > 0;
for k=1:dm(4)
    emu{k}       = emu{k}./semu;    
    emu{k}(~msk) = NaN;
end

% Compute derivatives
g = zeros([dm(1:2),3],'single');
H = zeros([dm(1:2),6],'single');
for k=1:dm(4)    
    % Gradient
    alpha = emu{k} - Z(:,:,k);
    for d=1:3
        g(:,:,d) = g(:,:,d) + alpha.*dmu{k,d};
    end

    % Hessian
    for k1=1:dm(4)
        if k1~=k
            tmp = -emu{k}.*emu{k1};
        else
            tmp = max(emu{k}.*(1-emu{k1}),0);
        end

        H(:,:,1)   = H(:,:,1) + tmp.*dmu{k,1}.*dmu{k1,1};
        H(:,:,2)   = H(:,:,2) + tmp.*dmu{k,2}.*dmu{k1,2};
        H(:,:,3)   = H(:,:,3) + tmp.*dmu{k,3}.*dmu{k1,3};
        H(:,:,4)   = H(:,:,4) + tmp.*dmu{k,1}.*dmu{k1,2};
        H(:,:,5)   = H(:,:,5) + tmp.*dmu{k,1}.*dmu{k1,3};
        H(:,:,6)   = H(:,:,6) + tmp.*dmu{k,2}.*dmu{k1,3};
    end
end

g(~isfinite(g)) = 0;
H(~isfinite(H)) = 0;
%==========================================================================