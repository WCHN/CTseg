function show_classes_in_Z(Z,dm,K)

Z    = reshape(Z,[dm K]);
step = ceil(0.05*dm(3));

figure(666)
for k=1:K
    subplot(1,K,k)
    img = [];
    for i=1:step:dm(3)
        R = Z(:,:,i,k);    
        img = [img; R];
    end
    imagesc(img); axis off xy
    drawnow
end
%==========================================================================