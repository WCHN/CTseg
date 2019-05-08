function [noise,mu_brain] = spm_noise_estimate_mod(Scans,speak,nr,nc,cnt_subplot)
% Estimate avarage noise from a series of images
% FORMAT noise = spm_noise_estimate(Scans)
% Scans - nifti structures or filenames of images
% noise - standard deviation estimate
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 2, speak       = 0; end
if nargin < 3, nr          = 0; end
if nargin < 4, nc          = 0; end
if nargin < 5, cnt_subplot = 0; end

if ~isa(Scans,'nifti'), Scans = nifti(Scans); end

noise    = zeros(numel(Scans),1);
mu_brain = zeros(numel(Scans),1);
for i=1:numel(Scans)
    Nii = Scans(i);
    f   = Nii.dat(:,:,:);
    if spm_type(Nii.dat.dtype(1:(end-3)),'intt')
        f(f==max(f(:))) = 0;
        x      = 0:Nii.dat.scl_slope:max(f(:));
        [h,x]  = hist(f(f>0),x);
    else
        x      = (0:1023)*(max(f(:))/1023);
        f(f==max(f(:))) = 0;
        [h,x]  = hist(f(f>0 & isfinite(f)),x);
    end
    [mg,nu,sd,mu] = spm_rice_mixture_mod(double(h(:)),double(x(:)),2,speak,nr,nc,cnt_subplot + i);
    
    noise(i)    = min(sd);
    
    x           = -nu.^2./(2*sd.^2);
    msk         = x>-20;
    Laguerre    = exp(x(msk)/2).*((1-x(msk)).*besseli(0,-x(msk)/2)-x(msk).*besseli(1,-x(msk)/2));
    Ey( msk)    = sqrt(pi*sd(msk).^2/2).*Laguerre;
    Ey(~msk)    = nu(~msk);
    mu_brain(i) = max(Ey);
end
%==========================================================================

%==========================================================================
function [mg,nu,sig,mu] = spm_rice_mixture_mod(h,x,K,speak,nr,nc,cnt_subplot)
% Fit a mixture of Ricians to a histogram
% FORMAT [mg,nu,sig] = rice_mixture(h,x,K)
% h    - histogram counts
% x    - bin positions (plot(x,h) to see the histogram)
% K    - number of Ricians
% mg   - integral under each Rician
% nu   - "mean" parameter of each Rician
% sig  - "standard deviation" parameter of each Rician
% mu   - "mean" parameter of each Rician, from sufficient statistics
%
% An EM algorithm is used, which involves alternating between computing
% belonging probabilities, and then the parameters of the Ricians.
% The Koay inversion technique is used to compute the Rician parameters
% from the sample means and standard deviations. This is described at
% http://en.wikipedia.org/wiki/Rician_distribution
%_______________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

mg  = ones(K,1)/K;
nu  = (0:(K-1))'*max(x)/(K+1);
sig = ones(K,1)*max(x)/(10*K);
mu  = zeros(K,1);

m0 = zeros(K,1);
m1 = zeros(K,1);
m2 = zeros(K,1);
ll = -Inf;
for iter=1:10000
    p  = zeros(numel(x),K);
    for k=1:K
        % Product Rule
        % p(class=k, x | mg, nu, sig) = p(class=k|mg) p(x | nu, sig, class=k)
        p(:,k) = mg(k)*ricepdf(x(:),nu(k),sig(k)^2);
    end

    % Sum Rule
    % p(x | mg, nu, sig) = \sum_k p(class=k, x | mg, nu, sig)
    sp  = sum(p,2)+eps;
    oll = ll;
    ll  = sum(log(sp).*h(:)); % Log-likelihood
    if ll-oll<1e-8*sum(h), break; end

%     fprintf('%g\n',ll);
%     md = mean(diff(x));
%     plot(x(:),p,'--',x(:),h/sum(h)/md,'b.',x(:),sp,'r'); drawnow

    % Bayes Rule
    % p(class=k | x, mg, nu, sig) = p(class=k, x | mg, nu, sig) / p(x | mg, nu, sig)
    p = bsxfun(@rdivide,p,sp);

    % Compute moments from the histograms, weighted by the responsibilities (p).
    for k=1:K
        m0(k) = sum(p(:,k).*h(:));              % Number of voxels in class k
        m1(k) = sum(p(:,k).*h(:).*x(:));        % Sum of the intensities in class k
        m2(k) = sum(p(:,k).*h(:).*x(:).*x(:));  % Sum of squares of intensities in class k
    end

    mg = m0/sum(m0); % Mixing proportions
    for k=1:K
        mu1 = m1(k)./m0(k);                                % Mean 
        mu2 = (m2(k)-m1(k)*m1(k)/m0(k)+1e-6)/(m0(k)+1e-6); % Variance

        % Compute nu & sig from mean and variance
        [nu(k),sig(k)] = moments2param(mu1,mu2);
        
        mu(k) = mu1;
    end
    %disp([nu'; sig'])
end

if speak
    if nr > 0
        subplot(nr,nc,cnt_subplot);
    end
    
    md = mean(diff(x));
    plot(x(:),p,'--',x(:),h/sum(h)/md,'b.',x(:),sp,'r'); drawnow
end

% % Background values p(x | z=bg)
% [~,ix] = min(sig);
% p1     = p(:,ix);
% 
% % Head values p(x | z=head)
% [~,ix] = max(sig);
% p2     = p(:,ix);
% 
% % Calculates intensity value for when p(head)=0.5 using
% % p(head) = p(x | z=head) / ( p(x | z=head) + p(x | z=bg) )
% pt = p2./(p2 + p1);
% 
% [~,ix]   = min(abs(pt - 0.5)); % Find closest value to 0.5
% val_brain = x(ix);
%==========================================================================

%==========================================================================
function [nu,sig] = moments2param(mu1,mu2)
% Rician parameter estimation (nu & sig) from mean (mu1) and variance
% (mu2) via the Koay inversion technique.
% This follows the scheme at
% https://en.wikipedia.org/wiki/Rice_distribution#Parameter_estimation_.28the_Koay_inversion_technique.29
% This Wikipedia description is based on:
% Koay, C.G. and Basser, P. J., Analytically exact correction scheme
% for signal extraction from noisy magnitude MR signals,
% Journal of Magnetic Resonance, Volume 179, Issue = 2, p. 317–322, (2006)

r     = mu1/sqrt(mu2);
theta = sqrt(pi/(4-pi));
if r>theta
    for i=1:256
        xi    = 2+theta^2-pi/8*exp(-theta^2/2)*((2+theta^2)*besseli(0,theta^2/4)+theta^2*besseli(1,theta^2/4))^2;
        g     = sqrt(xi*(1+r^2)-2);
        if abs(theta-g)<1e-6, break; end
        theta = g;
    end
    sig = sqrt(mu2)/sqrt(xi);
    nu  = sqrt(mu1^2+(xi-2)*sig^2);
else
    nu  = 0;
    sig = (2^(1/2)*(mu1^2 + mu2)^(1/2))/2;
end
%==========================================================================

%==========================================================================
function p = ricepdf(x,nu,sig2)
% Rician PDF
% p = ricepdf(x,nu,sig2)
% https://en.wikipedia.org/wiki/Rice_distribution#Characterization
p       = zeros(size(x));
tmp     = -(x.^2+nu.^2)./(2*sig2);
msk     = (tmp > -95) & (x*(nu/sig2) < 85) ; % Identify where Rice probability can be computed
p(msk)  = (x(msk)./sig2).*exp(tmp(msk)).*besseli(0,x(msk)*(nu/sig2)); % Use Rician distribution
p(~msk) = (1./sqrt(2*pi*sig2))*exp((-0.5/sig2)*(x(~msk)-nu).^2);      % Use Gaussian distribution
%==========================================================================