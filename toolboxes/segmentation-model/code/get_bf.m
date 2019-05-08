function [bf,ll] = get_bf(chan,dm,varargin)
% FORMAT [bf,ll] = get_bf(chan,dm,[bf,c,ll])
% chan - Structure array containing bias field encoding for all channels
% dm   - Dimension of the lattice on which to reconstruct the bias field
% c    - Channel to reconstruct (if not provided, reconstruct all channels)
% bf   - Reconstructed bias field
% ll   - Log-likelihood
%
% Reconstruct one or all bias fields on the full lattice from their DCT
% encoding.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
C  = numel(chan);
I  = prod(dm);
Iz = prod(dm(1:2));
nz = dm(3);   
if numel(varargin) == 0
    % Compute full bias field (for all channels)
    bf = zeros([I C],'single');
    ll = zeros(1,C);
    for c=1:C           
        ll(c) = double(-0.5*chan(c).T(:)'*chan(c).C*chan(c).T(:));    

        for z=1:nz
            ix = ix_slice(z,Iz);

            bf_c     = transf(chan(c).B1,chan(c).B2,chan(c).B3(z,:),chan(c).T);
            bf_c     = bf_c(:);      
            bf(ix,c) = single(exp(bf_c));        
        end
    end
else
    % Compute just for one channel
    bf = varargin{1};
    c  = varargin{2};
    ll = varargin{3};

    ll(c) = double(-0.5*chan(c).T(:)'*chan(c).C*chan(c).T(:));    

    for z=1:nz
        ix = ix_slice(z,Iz);

        bf_c     = transf(chan(c).B1,chan(c).B2,chan(c).B3(z,:),chan(c).T);
        bf_c     = bf_c(:);      
        bf(ix,c) = single(exp(bf_c));        
    end
end
%==========================================================================

%==========================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t  = zeros(size(B1,1),size(B2,1));
end
%==========================================================================