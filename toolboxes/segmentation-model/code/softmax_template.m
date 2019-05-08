function [mu,a1] = softmax_template(a,R,prop)
if nargin<3, prop = 0; end

d  = [size(a) 1 1 1];
mu = zeros([d(1:3),d(4)+1],'single');
if nargout==2
    a1 = zeros([d(1:3),d(4)+1],'single');
end

for j=1:size(a,3) % Loop over planes

    % Rotate the null-space back in to the data
    aj  = double(reshape(a(:,:,j,:),[d(1:2),d(4)]));
    sj  = zeros([d(1:2),d(4)+1]);
    for j1=1:d(4)+1
        sj(:,:,j1) = 0;
        for j2=1:d(4)
            sj(:,:,j1) = sj(:,:,j1) + R(j1,j2)*aj(:,:,j2);
        end
    end
    
    if nargout==2
        a1(:,:,j,:) = sj;
    end
    
    % Compute safe softmax
    sj          = bsxfun(@plus,sj,prop);
    mx_sj       = max(sj,[],3);
    sj          = exp(bsxfun(@minus,sj,mx_sj));
    s           = sum(sj,3);
    mu(:,:,j,:) = single(bsxfun(@rdivide,sj,s));
end
%==========================================================================