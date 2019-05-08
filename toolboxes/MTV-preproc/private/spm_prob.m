function varargout = spm_prob(varargin)
%__________________________________________________________________________
% Collection of tools for probability distributions (PDF, KL-div, ...).
%
% FORMAT out = spm_prob('Normal',  ...)
% FORMAT out = spm_prob('Gamma',   ...)
% FORMAT out = spm_prob('Wishart', ...)
% FORMAT out = spm_prob('Beta',    ...)
%
% FORMAT help spm_prob>function
% Returns the help file of the selected function.
%
% Equivalent keys
% ---------------
%   [n/normal/gaussian]
%   [g/gamma]
%   [ig/inverse-gamma]
%   [w/wishart]
%   [iw/inverse-wishart]
%   [ng/normal-gamma]
%   [nw/normal-wishart]
%   [ber/bernoulli]
%   [bin/binomial]
%   [c/categorical]
%   [m/multinomial]
%   [beta]
%   [d/dirichlet]
%--------------------------------------------------------------------------
% MISC
% ----
%
% FORMAT lg = spm_prob('LogGamma', a, p)
%   > Log of multivariate gamma function of order p
% FORMAT dg = spm_prob('DiGamma', a, p)
%   > Multivariate digamma function of order p
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

%--------------------------------------------------------------------------
% TODO
% ----
%
% - Functions returning only
%     * the parameter-dependent part
%     * the normalising constant
%   of the log-pdf.
% - Distributions:
%     * (Inv)-Wishart + Normal-(Inv)-Wishart
%     * Inv-Gamma     + Normal-(Inv)-Gamma
%     * Bernoulli/Categorical
%     * Dirichlet
%     * Laplace
% - Maximum-likelihood estimators (?)
% - Make functions work on arrays of observations
%   (for know it usually only works on scalar inputs)
%   It could allow using them for template update for exemple.
%   -> Needs standard inputs, especially in multivariate cases
%--------------------------------------------------------------------------
%
% If we want to extend these functions to volumes of observations or
% parameters, the input convention could be something along those lines:
%
% INPUT FORMAT
% ------------
% 
% - Single observations or single parameters should be scalar (or 1
%   dimensional vectors in the multivariate case).
% - Multiple observations should be of dimension > 3, with the 4th 
%   dimension (+5th in the multivariate case) being the feature space.
%   This allows us to deal with images or volumes in which each voxel is an
%   independent random variable.
% - Implicit expansion is performed if observations and parameters have
%   different dimensions.
% - Update functions need pre-computed ML estimators or sufficient
%   statistics. We do not deal with multiple observations of the same 
%   random variable.
%--------------------------------------------------------------------------

    if nargin == 0
        help spm_prob
        error('Not enough argument. Type ''help spm_prob'' for help.');
    end
    id = varargin{1};
    varargin = varargin(2:end);
    switch lower(id)
        case {'normal', 'n', 'gaussian'}
            [varargout{1:nargout}] = normal(varargin{:});
        case {'gamma', 'g'}
            [varargout{1:nargout}] = gamma(varargin{:});
        case {'wishart', 'w'}
            [varargout{1:nargout}] = wishart(varargin{:});
        case {'beta'}
            [varargout{1:nargout}] = beta(varargin{:});
        case {'loggamma'}
            [varargout{1:nargout}] = LogGamma(varargin{:});
        case {'digamma'}
            [varargout{1:nargout}] = DiGamma(varargin{:});
        otherwise
            help spm_prob
            error('Unknown function %s. Type ''help spm_prob'' for help.', id)
    end
end

%%
% =========================================================================
%   MISC
% =========================================================================

% -------------------------------------------------------------------------
function lg = LogGamma(a, p)
    if nargin < 2
        p = 1;
    end
    lg = (p*(p-1)/4)*log(pi);
    for i=1:p
        lg = lg + gammaln(a + (1-p)/2);
    end
end

% -------------------------------------------------------------------------

function dg = DiGamma(a, p, k)
    if nargin < 3
        k = 0;
        if nargin < 2
            p = 1;
        end
    end
    dg = 0;
    for i=1:p
        dg = dg + psi(k, a + (1-i)/2);
    end
end

%%
% =========================================================================
%   NORMAL
% =========================================================================

function varargout = normal(varargin)
%__________________________________________________________________________
% Characteristic functions of the (uni/multivariate) Normal distribution:
%   [pdf]       Probability density function
%   [ll/logpdf] Log-probability density function
%   [kl]        Kullback-Leibler divergence
%   [up/update] Conjugate (or Bayesian) update
%
%--------------------------------------------------------------------------
% General distribution
% --------------------
%
% The Normal distribution is parameterised by a mean parameter (mu) and
% a (co)variance (sigma) or precision (lambda) parameter. 
% For multivariate distributions, sigma is a KxK covariance matrix. 
% In the univariate case, it reduces to a scalar value, the variance.
%
% FORMAT pdf = spm_prob('Normal', 'pdf',    x, mu, sigma)
% FORMAT pdf = spm_prob('Normal', 'pdf',    x, mu, lambda, 'precision')
% FORMAT ll  = spm_prob('Normal', 'logpdf', x, mu, sigma)
% FORMAT ll  = spm_prob('Normal', 'logpdf', x, mu, lambda, 'precision')
%   >> (Log) Probability density function.
%
% FORMAT kl  = spm_prob('Normal', 'kl', mu1, sigma1,  mu0, sigma0)
% FORMAT kl  = spm_prob('Normal', 'kl', mu1, lambda1, mu0, lambda0, 'precision')
%   >> Kullback-Leibler divergence from N0 to N1 = KL(N1||N0)
%
%--------------------------------------------------------------------------
% Normal mean conjugate
% ---------------------
%
% The Normal distribution can be used as a conjugate prior for the mean
% parameter of another Normal distribution with known covariance.
% It is then parameterised by an expected mean (mu), a degrees of freedom 
% (n) and a known covariance (sigma) or precision (lambda). 
%
% FORMAT pdf = spm_prob('Normal', 'pdf',    x, mu, n, sigma)
% FORMAT pdf = spm_prob('Normal', 'pdf',    x, mu, n, lambda, 'precision')
% FORMAT ll  = spm_prob('Normal', 'logpdf', x, mu, n, sigma)
% FORMAT ll  = spm_prob('Normal', 'logpdf', x, mu, n, lambda, 'precision')
%   >> (Log) Probability density function.
%
% FORMAT kl  = spm_prob('Normal', 'kl', mu1, n1, mu0, n0, sigma)
% FORMAT kl  = spm_prob('Normal', 'kl', mu1, n1, mu0, n0, lambda, 'precision')
%   >> Kullback-Leibler divergence from N0 to N1 = KL(N1||N0)
%
% FORMAT [mu1, n1] = spm_prob('Normal', 'update', mu, n, mu0, n0)
%   >> Posterior parameters of the Normal distribution.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
    if nargin == 0
        help spm_prob>normal
        error('Not enough argument. Type ''help spm_prob>Normal'' for help.');
    end
    id = varargin{1};
    varargin = varargin(2:end);
    switch lower(id)
        case 'pdf'
            [varargout{1:nargout}] = normal_pdf(varargin{:});
        case {'logpdf', 'll'}
            [varargout{1:nargout}] = normal_logpdf(varargin{:});
        case 'kl'
            [varargout{1:nargout}] = normal_kl(varargin{:});
        case {'up', 'update'}
            [varargout{1:nargout}] = normal_up(varargin{:});
        case 'help'
            help spm_prob>normal
        otherwise
            help spm_prob>normal
            error('Unknown function %s. Type ''help spm_prob>Normal'' for help.', id)
    end
end

% -------------------------------------------------------------------------

function pdf = normal_pdf(x, mu, varargin)
% FORMAT pdf = spm_prob('Normal', 'pdf',    x, mu,    sigma)
% FORMAT pdf = spm_prob('Normal', 'pdf',    x, mu,    lambda, 'precision')
% FORMAT pdf = spm_prob('Normal', 'pdf',    x, mu, n, sigma)
% FORMAT pdf = spm_prob('Normal', 'pdf',    x, mu, n, lambda, 'precision')

    % Check if we are in the reparameterised case
    if nargin == 5
        if startsWith(varargin{3}, 'p', 'IgnoreCase', true);
            pdf = normal_pdf(x, mu, varargin{2}.*varargin{1}, 'precision');
        else
            pdf = normal_pdf(x, mu, varargin{2}./varargin{1});
        end
        return
    elseif nargin == 4 && ~ischar(varargin{2})
        pdf = normal_pdf(x, mu, varargin{2}./varargin{1});
        return
    end
    
    % Else, set default values
    if nargin < 4
        mode = 'covariance';
        if nargin < 2
            mu = zeros(size(x));
        end
        if nargin < 3
            varargin{1} = eye(numel(mu));
        end
    else
        mode = varargin{2};
    end
    precision = startsWith(mode, 'p', 'IgnoreCase', true);
    
    % Usual PDF
    if precision
        pdf = (exp(spm_matcomp('LogDet', varargin{1}/(2*pi))))^(0.5) * exp(-0.5*(x(:)-mu(:))'*varargin{1}*(x(:)-mu(:)));
    else
        pdf = (exp(spm_matcomp('LogDet', varargin{1}*2*pi)))^(-0.5) * exp(-0.5*(x(:)-mu(:))'*(varargin{1}\(x(:)-mu(:))));
    end
    
end

% -------------------------------------------------------------------------

function pdf = normal_logpdf(x, mu, varargin)
% FORMAT ll  = spm_prob('Normal', 'logpdf', x, mu,    sigma)
% FORMAT ll  = spm_prob('Normal', 'logpdf', x, mu,    lambda, 'precision')
% FORMAT ll  = spm_prob('Normal', 'logpdf', x, mu, n, sigma)
% FORMAT ll  = spm_prob('Normal', 'logpdf', x, mu, n, lambda, 'precision')

    % Check if we are in the reparameterised case
    if nargin == 5
        if startsWith(varargin{3}, 'p', 'IgnoreCase', true);
            pdf = normal_pdf(x, mu, varargin{2}*varargin{1}, 'precision');
        else
            pdf = normal_pdf(x, mu, varargin{2}/varargin{1});
        end
        return
    elseif nargin == 4 && ~ischar(varargin{2})
        pdf = normal_pdf(x, mu, varargin{2}/varargin{1});
        return
    end
    
    % Else, set default values
    if nargin < 4
        mode = 'covariance';
        if nargin < 2
            mu = zeros(size(x));
        end
        if nargin < 3
            varargin{1} = eye(numel(mu));
        end
    else
        mode = varargin{2};
    end
    precision = startsWith(mode, 'p', 'IgnoreCase', true);
    K = size(varargin{1}, 1);
    
    if precision
        pdf = -0.5*( K*log(2*pi) - spm_matcomp('LogDet', varargin{1}) + (x(:)-mu(:))'*varargin{1}*(x(:)-mu(:)) );
    else
        pdf = -0.5*( K*log(2*pi) + spm_matcomp('LogDet', varargin{1}) + (x(:)-mu(:))'*(varargin{1}\(x(:)-mu(:))) );
    end
end

% -------------------------------------------------------------------------

function kl = normal_kl(mu1, par1, mu0, par0, varargin)
% FORMAT kl = spm_prob('Normal', 'kl', mu1, sigma1,  mu0, sigma0)
% FORMAT kl = spm_prob('Normal', 'kl', mu1, lambda1, mu0, lambda0, 'precision')
% FORMAT kl = spm_prob('Normal', 'kl', mu1, n1,      mu0, n0,      sigma)
% FORMAT kl = spm_prob('Normal', 'kl', mu1, n1,      mu0, n0,      lambda, 'precision')

    % Check if we are in the reparameterised case
    if nargin == 6
        K = size(varargin{1}, 1);
        if startsWith(varargin{2}, 'p', 'IgnoreCase', true);
            kl = 0.5 * ( K * (par0/par1 - log(par0/par1) - 1) ...
                         + par0 * (mu0(:)-mu1(:))'*varargin{1}*(mu0(:)-mu1(:)) );
        else
            kl = 0.5 * ( K * (par0/par1 - log(par0/par1) - 1) ...
                         + par0 * (mu0(:)-mu1(:))'*(varargin{1}\(mu0(:)-mu1(:))) );
        end
        return
    elseif nargin == 5 && ~ischar(varargin{1})
        K  = size(varargin{1}, 1);
        kl = 0.5 * ( K * (par0/par1 - log(par0/par1) - 1) ...
                     + par0 * (mu0(:)-mu1(:))'*(varargin{1}\(mu0(:)-mu1(:))) );
        return
    end
    
    % Else, set default values
    if nargin < 5
        mode = 'covariance';
    else
        mode = varargin{1};
    end
    precision = startsWith(mode, 'p', 'IgnoreCase', true);
    K = size(par1, 1);
    
    % Common KL-divergence
    if precision
        kl = 0.5 * ( trace(par1\par0) ...
                     - spm_matcomp('LogDet', par0) ...
                     + spm_matcomp('LogDet', par1) ...
                     - K ...
                     + (mu0(:)-mu1(:))'*par0*(mu0(:)-mu1(:)) );
    else
        kl = 0.5 * ( trace(par0\par1) ...
                     - spm_matcomp('LogDet', par1) ...
                     + spm_matcomp('LogDet', par0) ...
                     - K ...
                     + (mu0(:)-mu1(:))'*(par1\(mu0(:)-mu1(:))) );
    end
    
end

% -------------------------------------------------------------------------

function [mu1, n1] = normal_up(mu, n, mu0, n0)
% FORMAT [mu1, n1] = spm_prob('Normal', 'update', mu, n, mu0, n0)

    n1 = n + n0;
    mu1 = (n0 * mu0 + n * mu)/n1;

end

%%
% =========================================================================
%   GAMMA
% =========================================================================

function varargout = gamma(varargin)
%__________________________________________________________________________
% Characteristic functions of the Gamma distribution:
%   [pdf]       Probability density function
%   [ll/logpdf] Log-probability density function
%   [kl]        Kullback-Leibler divergence
%   [up/update] Conjugate (or Bayesian) update
%   [E]         Expected value (E[x])
%   [Elog]      Expected log (E[ln x])
%   [V]         Variance (V[x])
%   [Vlog]      Variance of the log (V[ln x])
%
% The Gamma distribution is a conjugate prior for a Normal precision (or
% precision magnitude) with known mean, for a Gamma rate with known shape 
% or in general for any rate parameter of an Exponential family
% distribution.
%
%--------------------------------------------------------------------------
% General distribution
% --------------------
%
% The Gamma distribution is parameterised by a shape parameter (alpha) and
% a rate parameter (beta).
%
% FORMAT pdf = spm_prob('Gamma', 'pdf',    x, alpha, beta)
% FORMAT ll  = spm_prob('Gamma', 'logpdf', x, alpha, alpha)
%   >> (Log) Probability density function.
%
% FORMAT e  = spm_prob('Gamma', 'E',    alpha, beta)
% FORMAT el = spm_prob('Gamma', 'Elog', alpha, beta)
% FORMAT v  = spm_prob('Gamma', 'V',    alpha, beta)
% FORMAT vl = spm_prob('Gamma', 'Vlog', alpha, beta)
%   >> Mean and variance
%
% FORMAT kl  = spm_prob('Gamma', 'kl', alpha1, beta1, alpha0, beta0)
%   >> Kullback-Leibler divergence from G0 to G1 = KL(G1||G0).
%
%--------------------------------------------------------------------------
% Normal precision conjugate
% --------------------------
%
% The Gamma distribution is parameterised by a mean precision parameter 
% (lambda) and a degrees of freedom (n). It can be a precision *magnitude*
% parameter of K > 1.
%
% FORMAT pdf = spm_prob('Gamma', 'pdf',    x, lambda, n, K, ('normal'))
% FORMAT ll  = spm_prob('Gamma', 'logpdf', x, lambda, n, K, ('normal'))
%   >> (Log) Probability density function.
%
% FORMAT e  = spm_prob('Gamma', 'E',    lambda, n, K, ('normal'))
% FORMAT el = spm_prob('Gamma', 'Elog', lambda, n, K, ('normal'))
% FORMAT v  = spm_prob('Gamma', 'V',    lambda, n, K, ('normal'))
% FORMAT vl = spm_prob('Gamma', 'Vlog', lambda, n, K, ('normal'))
%   >> Mean and variance.
%
% FORMAT kl  = spm_prob('Gamma', 'kl', lam1, n1, lam0, n0, K, ('normal'))
%   >> Kullback-Leibler divergence from G0 to G1 = KL(G1||G0).
%
% FORMAT [lam1, n1] = spm_prob('Gamma', 'up', lam, n,     lam0, n0, ('normal'))
% FORMAT [lam1, n1] = spm_prob('Gamma', 'up', s0, s1, s2, lam0, n0, 
%                                             (mu=0), (Lam=eye), ('normal'))
%   >> Posterior parameters of the Gamma distribution.
%
%--------------------------------------------------------------------------
% Gamma rate conjugate
% --------------------
%
% The Gamma distribution is parameterised by a mean rate parameter 
% (beta), a degrees of freedom (n), and a known shape parameter (alpha).
%
% FORMAT pdf = spm_prob('Gamma', 'pdf',    x, beta, n, alpha, 'gamma')
% FORMAT ll  = spm_prob('Gamma', 'logpdf', x, beta, n, alpha, 'gamma')
%   >> (Log) Probability density function.
%
% FORMAT e  = spm_prob('Gamma', 'E',    beta, n, alpha, 'gamma')
% FORMAT el = spm_prob('Gamma', 'Elog', beta, n, alpha, 'gamma')
% FORMAT v  = spm_prob('Gamma', 'V',    beta, n, alpha, 'gamma')
% FORMAT vl = spm_prob('Gamma', 'Vlog', beta, n, alpha, 'gamma')
%   >> Mean and variance.
%
% FORMAT kl  = spm_prob('Gamma', 'kl', beta1, n1, beta0, n0, alpha, 'gamma')
%   >> Kullback-Leibler divergence from G0 to G1 = KL(G1||G0).
%
% FORMAT [beta1, n1] = spm_prob('Gamma', 'up', beta, n, beta0, n0, 'gamma')
% FORMAT [beta1, n1] = spm_prob('Gamma', 'up', s0, s1,  beta0, n0, alpha, 'gamma')
%   >> Posterior parameters of the Gamma distribution.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
    if nargin == 0
        help spm_prob>gamma
        error('Not enough argument. Type ''help spm_prob>Gamma'' for help.');
    end
    id = varargin{1};
    varargin = varargin(2:end);
    switch lower(id)
        case 'pdf'
            [varargout{1:nargout}] = gamma_pdf(varargin{:});
        case {'logpdf', 'll'}
            [varargout{1:nargout}] = gamma_logpdf(varargin{:});
        case 'kl'
            [varargout{1:nargout}] = gamma_kl(varargin{:});
        case {'up', 'update'}
            [varargout{1:nargout}] = gamma_up(varargin{:});
        case 'e'
            [varargout{1:nargout}] = gamma_e(varargin{:});
        case 'v'
            [varargout{1:nargout}] = gamma_v(varargin{:});
        case 'elog'
            [varargout{1:nargout}] = gamma_elog(varargin{:});
        case 'vlog'
            [varargout{1:nargout}] = gamma_vlog(varargin{:});
        case 'help'
            help spm_prob>gamma
        otherwise
            help spm_prob>gamma
            error('Unknown function %s. Type ''help spm_prob>Gamma'' for help.', id)
    end
end

% -------------------------------------------------------------------------

function pdf = gamma_pdf(x, varargin)
% FORMAT pdf = gamma_pdf(x, alpha,  beta)
% FORMAT pdf = gamma_pdf(x, lambda, n,    K,     ('normal'))
% FORMAT pdf = gamma_pdf(x, beta,   n,    alpha, 'gamma')

    % Check if we are in the reparameterised case
    if nargin == 4
        if ischar(varargin{3})
            varargin{3} = 1;
        end
        a   = 0.5 * bsxfun(@times, varargin{3}, varargin{2});
        b   = bsxfun(@rdivide, a, varargin{1});
    elseif nargin > 4 
        if startsWith(varargin{4}, 'n', 'IgnoreCase', true)
            a   = 0.5 * bsxfun(@times, varargin{3}, varargin{2});
            b   = bsxfun(@rdivide, a, varargin{1});
        elseif startsWith(varargin{4}, 'g', 'IgnoreCase', true)
            a   = bsxfun(@times, varargin{3}, varargin{2});
            b   = bsxfun(@rdivide, a, varargin{1});
        end
    else
        a = varargin{1};
        b = varargin{2};
    end
    
    % Usual pdf
    pdf = gampdf(x, a, 1./b);
    
end

% -------------------------------------------------------------------------

function pdf = gamma_logpdf(x, varargin)
% FORMAT pdf = gamma_logpdf(x, alpha,  beta)
% FORMAT pdf = gamma_logpdf(x, lambda, n,    K,     ('normal'))
% FORMAT pdf = gamma_logpdf(x, beta,   n,    alpha, 'gamma')

    % Check if we are in the reparameterised case
    if nargin == 4
        if ischar(varargin{3})
            varargin{3} = 1;
        end
        a   = 0.5 * bsxfun(@times, varargin{3}, varargin{2});
        b   = bsxfun(@rdivide, a, varargin{1});
    elseif nargin > 4 
        if startsWith(varargin{4}, 'n', 'IgnoreCase', true)
            a   = 0.5 * bsxfun(@times, varargin{3}, varargin{2});
            b   = bsxfun(@rdivide, a, varargin{1});
        elseif startsWith(varargin{4}, 'g', 'IgnoreCase', true)
            a   = bsxfun(@times, varargin{3}, varargin{2});
            b   = bsxfun(@rdivide, a, varargin{1});
        end
    else
        a = varargin{1};
        b = varargin{2};
    end
    
    % Usual log-pdf
    pdf = bsxfun(@plus, bsxfun(@times, a, log(b)), ...
                        bsxfun(@times, a-1, log(x)));
    pdf = bsxfun(@minus, pdf, bsxfun(@times, b, x));
    pdf = bsxfun(@minus, pdf, gammaln(a));
    
end

% -------------------------------------------------------------------------

function kl = gamma_kl(varargin)
% FORMAT pdf = gamma_kl(alpha1,  beta1, alpha0,  beta0)
% FORMAT pdf = gamma_kl(lambda1, n1,    lambda0, n0,    K,     ('normal'))
% FORMAT pdf = gamma_kl(beta1,   n1,    beta0,   n0,    alpha, 'gamma')

    % Check if we are in the reparameterised case
    if nargin == 5
        if ischar(varargin{5})
            varargin{5} = 1;
        end
        a1  = 0.5*bsxfun(@times, varargin{5}, varargin{2});
        b1  = bsxfun(@rdivide, a1, varargin{1});
        a0  = 0.5*bsxfun(@times, varargin{5}, varargin{4});
        b0  = bsxfun(@rdivide, a0, varargin{3});
    elseif nargin > 5
        if startsWith(varargin{6}, 'n', 'IgnoreCase', true)
            a1  = 0.5*bsxfun(@times, varargin{5}, varargin{2});
            b1  = bsxfun(@rdivide, a1, varargin{1});
            a0  = 0.5*bsxfun(@times, varargin{5}, varargin{4});
            b0  = bsxfun(@rdivide, a0, varargin{3});
        elseif startsWith(varargin{6}, 'g', 'IgnoreCase', true)
            a1  = bsxfun(@times, varargin{5}, varargin{2});
            b1  = bsxfun(@rdivide, a1, varargin{1});
            a0  = bsxfun(@times, varargin{5}, varargin{4});
            b0  = bsxfun(@rdivide, a0, varargin{3});
        end
    else
        a1 = varargin{1};
        b1 = varargin{2};
        a0 = varargin{3};
        b0 = varargin{4};  
    end
    
    % Usual KL
    % KL = - a0*log(b0/b1) + a1*(b0/b1 - 1) + (a1-a0)*psi(a1)
    %      + ln G(a0) - ln G(a1)
    kl = bsxfun(@times, bsxfun(@minus, a1, a0), psi(a1));
    kl = bsxfun(@minus, kl, bsxfun(@times, a0, log(bsxfun(@rdivide, b0, b1))));
    kl = bsxfun(@minus, kl, bsxfun(@times, a1, bsxfun(@rdivide, b0, b1)-1));
    kl = bsxfun(@plus,  kl, gammaln(a0)) ;
    kl = bsxfun(@minus, kl, gammaln(a1)) ;
    
end

% -------------------------------------------------------------------------

function [par1, n1] = gamma_up(varargin)
% FORMAT [lam1, n1]  = gamma_up(lam, n,        lam0, n0, ('normal'))
% FORMAT [lam1, n1]  = gamma_up(ss0, ss1, ss2, lam0, n0, (mu=0), (Lambda=eye), 'normal')
% FORMAT [beta1, n1] = gamma_up(beta, n, beta0, n0, 'gamma')
% FORMAT [beta1, n1] = gamma_up(ss0, ss1,  beta0, n0, alpha, 'gamma')

    if ischar(varargin{end}) && ...
       startsWith(varargin{end}, 'g', 'IgnoreCase', true)
        % -----
        % GAMMA
        % -----
        varargin = varargin(1:end-1);
        
        if nargin > 4
            % Sufficient statistics case
            ss0   = varargin{1};
            ss1   = varargin{2};
            beta0 = varargin{3};
            n0    = varargin{4};
            alpha = varargin{5};
            n1 = n0 + ss0;
            par1 = n0/beta0 + ss1/alpha;
            par1 = n1 / par1;
        else
            % Average case
            beta    = varargin{1};
            n       = varargin{2};
            beta0   = varargin{3};
            n0      = varargin{4};
            n1      = n + n0;
            par1 = n1 / (n0/beta0 + n/beta);
        end
        
    else
        % ------
        % NORMAL
        % ------
        if ischar(varargin{end})
            varargin = varargin(1:end-1);
        end
    
        if nargin > 4
            % Sufficient statistics case
            ss0     = varargin{1};
            ss1     = varargin{2};
            ss2     = varargin{3};
            lambda0 = varargin{4};
            n0      = varargin{5};
            K  = size(ss2, 1);
            if nargin < 7
                Lambda = eye(K);
                if nargin < 6
                    mu = zeros(size(ss1));
                else
                    mu = varargin{6};
                end
            else
                Lambda = varargin{7};
            end
            n      = ss0;
            lambda = trace(ss2*Lambda) ...
                     - 2 * mu'*Lambda*ss1 ...
                     + ss0 * mu'*Lambda*mu;
            lambda = K*ss0/lambda;
            [par1, n1] = gamma_up(lambda, n, lambda0, n0);
        else
            % Average case
            lambda  = varargin{1};
            n       = varargin{2};
            lambda0 = varargin{3};
            n0      = varargin{4};
            n1      = n + n0;
            par1 = n1 / (n0/lambda0 + n/lambda);
        end
    
    end
    
end

% -------------------------------------------------------------------------

function out = gamma_e(varargin)
% FORMAT e = gamma_e(a, b)
% FORMAT e = gamma_e(lam, n, (K), ('normal'))
% FORMAT e = gamma_e(b,   n, a,    'gamma')

    if ischar(varargin{end}) && ...
       startsWith(varargin{end}, 'g', 'IgnoreCase', true)
        % ----------
        % GAMMA CONJ
        % ----------
            out = varargin{2};
            
    elseif ( ischar(varargin{end}) && ...
             startsWith(varargin{end}, 'n', 'IgnoreCase', true) ) || ...
           (nargin == 3)
        % -----------
        % NORMAL CONJ
        % -----------
            out = varargin{2};
            
    else
        % -----
        % GAMMA
        % -----
            out = bsxfun(@rdivide, varargin{1}, varargin{2});
    end
end

% -------------------------------------------------------------------------

function out = gamma_v(varargin)
% FORMAT v = gamma_v(a, b)
% FORMAT v = gamma_v(lam, n, (K), ('normal'))
% FORMAT v = gamma_v(b,   n, a,    'gamma')

    if ischar(varargin{end}) && ...
       startsWith(varargin{end}, 'g', 'IgnoreCase', true)
        % ----------
        % GAMMA CONJ
        % ----------
        a = bsxfun(@times, varargin{3}, varargin{2});
        b = bsxfun(@rdivide, a, varargin{1});
                             
    elseif ( ischar(varargin{end}) && ...
             startsWith(varargin{end}, 'n', 'IgnoreCase', true) ) || ...
           (nargin == 3)
        % -----------
        % NORMAL CONJ
        % -----------
        if ischar(varargin{3})
            varargin{3} = 1;
        end
        a = 0.5 * bsxfun(@times, varargin{3}, varargin{2});
        b = bsxfun(@rdivide, a, varargin{1});
                                 
    else
        % -----
        % GAMMA
        % -----
        a = varargin{1};
        b = varargin{2};
        
    end
    out = bsxfun(@rdivide, a, b.^2);
end

% -------------------------------------------------------------------------

function out = gamma_elog(varargin)
% FORMAT el = gamma_elog(a, b)
% FORMAT el = gamma_elog(lam, n, (K), ('normal'))
% FORMAT el = gamma_elog(b,   n, a,    'gamma')

    if ischar(varargin{end}) && ...
       startsWith(varargin{end}, 'g', 'IgnoreCase', true)
        % ----------
        % GAMMA CONJ
        % ----------
        a = bsxfun(@times, varargin{3}, varargin{2});
        b = bsxfun(@rdivide, a, varargin{1});
            
    elseif ( ischar(varargin{end}) && ...
             startsWith(varargin{end}, 'n', 'IgnoreCase', true) ) || ...
           (nargin == 3)
        % -----------
        % NORMAL CONJ
        % -----------
        if ischar(varargin{3})
            varargin{3} = 1;
        end
        a = 0.5 * bsxfun(@times, varargin{3}, varargin{2});
        b = bsxfun(@rdivide, a, varargin{1});
              
    else
        % -----
        % GAMMA
        % -----
        a = varargin{1};
        b = varargin{2};
    end
    
    out = bsxfun(@minus, psi(a), log(b));
end

% -------------------------------------------------------------------------

function out = gamma_vlog(varargin)
% FORMAT vl = gamma_vlog(a, b)
% FORMAT vl = gamma_vlog(lam, n, (K), ('normal'))
% FORMAT vl = gamma_vlog(b,   n, a,    'gamma')

    if ischar(varargin{end}) && ...
       startsWith(varargin{end}, 'g', 'IgnoreCase', true)
        % ----------
        % GAMMA CONJ
        % ----------
        a = bsxfun(@times, varargin{3}, varargin{2});
        
    elseif ( ischar(varargin{end}) && ...
             startsWith(varargin{end}, 'n', 'IgnoreCase', true) ) || ...
           (nargin == 3)
        % -----------
        % NORMAL CONJ
        % -----------
        if ischar(varargin{3})
            varargin{3} = 1;
        end
        a = 0.5 * bsxfun(@times, varargin{3}, varargin{2});
        
    else
        % -----
        % GAMMA
        % -----
        a = varargin{1};
    end
    
    out = psi(1,a);
end

%%
% =========================================================================
%   WISHART
% =========================================================================

function varargout = wishart(varargin)
%__________________________________________________________________________
% Characteristic functions of the Wishart distribution:
%   [pdf]       Probability density function
%   [ll/logpdf] Log-probability density function
%   [kl]        Kullback-Leibler divergence
%   [up/update] Conjugate (or Bayesian) update
%   [E]         Expected value (E[X])
%   [Elogdet]   Expected log determinant (E[ln det X])
%   [V]         Variance (V[X])
%   [Vlogdet]   Variance of the log determinant (V[ln det X])
%
% The Wishart distribution is a conjugate prior for a multivariate Normal 
% precision matrix with known mean.
%
%--------------------------------------------------------------------------
% General distribution
% --------------------
%
% The Wishart distribution is parameterised by a scale matrix (V) and a
% degrees of freedom (n). It can be seen as the distribution of the sum of
% n independent centered multivariate Normal variables with precision 
% matrix V.
%
% FORMAT pdf = spm_prob('Wishart', 'pdf',    X, V, n)
% FORMAT ll  = spm_prob('Wishart', 'logpdf', X, V, n)
%   >> (Log) Probability density function.
%
% FORMAT e    = spm_prob('Wishart', 'E',       V, n)
% FORMAT el   = spm_prob('Wishart', 'Elogdet', V, n)
% FORMAT v    = spm_prob('Wishart', 'V',       V, n)
% FORMAT vl   = spm_prob('Wishart', 'Vlogdet', V, n)
% FORMAT logZ = spm_prob('Wishart', 'logZ',    V, n)
%   >> Mean and variance
%
% FORMAT kl  = spm_prob('Wishart', 'kl', V1, n1, V0, n0)
%   >> Kullback-Leibler divergence from W0 to W1 = KL(W1||W0).
%
%--------------------------------------------------------------------------
% Normal precision matrix conjugate
% ---------------------------------
%
% The Wishart distribution is parameterised by a mean precision parameter 
% (Lambda) and a degrees of freedom (n).
%
% FORMAT pdf = spm_prob('Wishart', 'pdf',    X, Lambda, n, 'normal')
% FORMAT ll  = spm_prob('Wishart', 'logpdf', X, Lambda, n, 'normal')
%   >> (Log) Probability density function.
%
% FORMAT e  = spm_prob('Wishart', 'E',       Lambda, n, 'normal')
% FORMAT el = spm_prob('Wishart', 'Elogdet', Lambda, n, 'normal')
% FORMAT v  = spm_prob('Wishart', 'V',       Lambda, n, 'normal')
% FORMAT vl = spm_prob('Wishart', 'Vlogdet', Lambda, n, 'normal')
%   >> Mean and variance.
%
% FORMAT kl  = spm_prob('Wishart', 'kl', Lam1, n1, Lam0, n0, 'normal')
%   >> Kullback-Leibler divergence from W0 to W1 = KL(W1||W0).
%
% FORMAT [lam1, n1] = spm_prob('Wishart', 'up', Lam, n,     Lam0, n0)
% FORMAT [lam1, n1] = spm_prob('Wishart', 'up', s0, s1, s2, Lam0, n0, (mu=0))
%   >> Posterior parameters of the Wishart distribution.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
    if nargin == 0
        help spm_prob>wishart
        error('Not enough argument. Type ''help spm_prob>Wishart'' for help.');
    end
    id = varargin{1};
    varargin = varargin(2:end);
    switch lower(id)
        case 'pdf'
            [varargout{1:nargout}] = wishart_pdf(varargin{:});
        case {'logpdf', 'll'}
            [varargout{1:nargout}] = wishart_logpdf(varargin{:});
        case 'kl'
            [varargout{1:nargout}] = wishart_kl(varargin{:});
        case {'up', 'update'}
            [varargout{1:nargout}] = wishart_up(varargin{:});
        case 'e'
            [varargout{1:nargout}] = wishart_e(varargin{:});
        case 'elogdet'
            [varargout{1:nargout}] = wishart_elogdet(varargin{:});
        case 'vloget'
            [varargout{1:nargout}] = wishart_vlogdet(varargin{:});
        case 'logz'
            [varargout{1:nargout}] = wishart_logZ(varargin{:});            
        case 'help'
            help spm_prob>wishart
        otherwise
            help spm_prob>wishart
            error('Unknown function %s. Type ''help spm_prob>Wishart'' for help.', id)
    end
end

% -------------------------------------------------------------------------

function pdf = wishart_pdf(X, varargin)
% FORMAT pdf = wishart_pdf(X, V,      n)
% FORMAT pdf = wishart_pdf(X, Lambda, n, 'normal')

    % Check if we are in the reparameterised case
    if nargin == 4
        pdf = wishart_pdf(X, varargin{1}/varargin{2}, varargin{2});
        return
    end
    
    % Usual pdf
    V   = varargin{1};
    n   = varargin{2};
    K   = size(V, 1);
    pdf = det(X)^(n-K-1) * exp(-0.5*trace(V\X)) ...
          / ( 2^(n*K/2) * det(V)^(n/2) * exp(LogGamma(0.5*n, K)) );
    
end

% -------------------------------------------------------------------------

function pdf = wishart_logpdf(X, varargin)
% FORMAT pdf = wishart_logpdf(X, V,      n)
% FORMAT pdf = wishart_logpdf(X, Lambda, n, 'normal')

    % Check if we are in the reparameterised case
    if nargin == 4
        pdf = wishart_logpdf(X, varargin{1}/varargin{2}, varargin{2});
        return
    end
    
    % Usual pdf
    V   = varargin{1};
    n   = varargin{2};
    K   = size(V, 1);
    pdf =   0.5*(n-K-1)*spm_matcomp('LogDet', X) ...
          - 0.5*trace(V\X) ...
          - 0.5*n*K*log(2) ...
          - 0.5*n*spm_matcomp('LogDet', V) ...
          - LogGamma(0.5*n, K);
    
end

% -------------------------------------------------------------------------

function kl = wishart_kl(varargin)
% FORMAT kl = wishart_kl(V1,      n1, V0,      n0)
% FORMAT kl = wishart_kl(lambda1, n1, lambda0, n0, 'normal')

    % Check if we are in the reparameterised case
    if nargin == 5
        kl = wishart_kl(varargin{1}/varargin{2}, varargin{2}, ...
                      varargin{3}/varargin{4}, varargin{4});
        return
    end
    
    % Usual KL
    V1 = varargin{1};
    n1 = varargin{2};
    V0 = varargin{3};
    n0 = varargin{4};
    K  = size(V1, 1);
    kl =   0.5*n0*(spm_matcomp('LogDet', V0) - spm_matcomp('LogDet', V1)) ...
         + 0.5*n1*(trace(V0\V1) - K) ...
         + 0.5*(n1 - n0)*DiGamma(0.5*n1, K) ...
         + LogGamma(0.5*n0, K) - LogGamma(0.5*n1, K);
    
end

% -------------------------------------------------------------------------

function [Lam1, n1] = wishart_up(varargin)
% FORMAT [Lam1, n1]  = wishart_up(Lam, n,        Lam0, n0)
% FORMAT [Lam1, n1]  = wishart_up(ss0, ss1, ss2, Lam0, n0, (mu=0))

    
    if nargin > 4
        % Sufficient statistics case
        ss0  = varargin{1};
        ss1  = varargin{2};
        ss2  = varargin{3};
        Lam0 = varargin{4};
        n0   = varargin{5};
        if nargin < 6
            iV = ss2;
        else
            iV = ss2 - 2*ss1*mu' + ss0*mu*mu';
        end
        n1   = n0+ss0;
        if n0, Lam1 = (n0*inv(Lam0) + iV)/n1;
        else,  Lam1 = iV/n1; end
        % stable inverse
        Lam1 = spm_matcomp('Inv', Lam1);
    else
        % Average case
        Lam  = varargin{1};
        n    = varargin{2};
        Lam0 = varargin{3};
        n0   = varargin{4};
        n1   = n + n0;
        if n0, Lam1 = n1*spm_matcomp('Inv', ...
                            n0 * spm_matcomp('Inv',Lam0) + ...
                            n  * spm_matcomp('Inv',Lam));
        else,  Lam1 = Lam; end
    end
    
    
end

% -------------------------------------------------------------------------

function out = wishart_e(V, n, mode)
    if nargin < 3 || mode(1) ~= 'n'
        out = n*V;
    else
        out = V;
    end
end

% -------------------------------------------------------------------------

function out = wishart_elogdet(V, n, mode)
    if nargin < 3 || mode(1) ~= 'n'
        K   = size(V, 1);
        out = DiGamma(0.5*n, K) + K*log(2) + spm_matcomp('LogDet', V);
    else
        out = wishart_elogdet(V/n, n);
    end
end

% -------------------------------------------------------------------------

function out = wishart_vlogdet(V, n, ~)
    K = size(V, 1);
    out = 0;
    for i=1:K
        out = out + psi(1, 0.5*(n+1-i));
    end
end

% -------------------------------------------------------------------------

function out = wishart_logZ(V,n)
    M     = size(V,1);
    if M >= n+1
       %warning('SPM:Wishart','Can not normalise a Wishart distribution (M=%d, nu=%f)', M,nu);
        out = 0;
        return;
    end
    lGamM = M*(M-1)/4*log(pi);
    for m=1:M
        lGamM = lGamM + gammaln((n + 1 - m)/2);
    end
    out = 0.5*n*(spm_matcomp('LogDet',V) + M*log(2)) + lGamM;
end

%%
% =========================================================================
%   BETA
% =========================================================================

function varargout = beta(varargin)
%__________________________________________________________________________
% Characteristic functions of the Beta distribution:
%
%   [pdf]       Probability density function
%   [ll/logpdf] Log-probability density function
%   [kl]        Kullback-Leibler divergence
%   [up/update] Conjugate (or Bayesian) update
%   [E]         Expected value (E[x])
%   [Elog]      Expected log value (E[ln x])
%   [V]         Variance (V[x])
%   [Vlog]      Variance of the log value (V[ln x])
%
%--------------------------------------------------------------------------
% Relationship to other distributions
% -----------------------------------
%
% The Beta distribution is a conjugate prior for the proportion parameter 
% of the Bernoulli, Binomial, negative Binomial and Geometric 
% distributions. 
%
% The beta prameter is a special case of the Dirichlet distribution when
% the number of classes is K = 2, i.e., the underlying trial is binary.
%
%--------------------------------------------------------------------------
% General distribution
% --------------------
%
% The Beta distribution is parameterised by two strictly positive shape 
% parameters, a and b, and is defined over [0,1]. Parameters a and b can be
% seen as concentration parameters, i.e., they reflect both class
% proportions and their variance.
% The expected value of the distribution is E[x] = a/(a+b)
% The variance of the distribution is       V[x] = ab/[(a+b)^2*(a+b+1)]
%
% FORMAT pdf = spm_prob('Beta', 'pdf',    x, a, b)
% FORMAT ll  = spm_prob('Beta', 'logpdf', x, a, b)
%   >> (Log) Probability density function.
%
% FORMAT e  = spm_prob('Beta', 'E',    a, b)
% FORMAT el = spm_prob('Beta', 'Elog', a, b)
% FORMAT v  = spm_prob('Beta', 'V',    a, b)
% FORMAT vl = spm_prob('Beta', 'Vlog', a, b)
%   >> Mean and variance
%
% FORMAT kl  = spm_prob('Beta', 'kl', a1, b1, a0, b0)
%   >> Kullback-Leibler divergence from B0 to B1 = KL(B1||B0).
%
%--------------------------------------------------------------------------
% Bernoulli probability conjugate
% -------------------------------
%
% To make the distribution easier to manipulate when used as a conjugate
% prior for a Bernoulli distribution, we reparameterise it by setting:
% - a = np
% - b = n(1-p)
% The expected value of the distribution becomes E[x] = p
% The variance of the distribution becomes       V[x] = p(1-p)/(n+1)
% Then, the Beta distribution is parameterised by a mean probability  
% parameter (p) and a degrees of freedom (n).
%
% FORMAT pdf = spm_prob('Beta', 'pdf',    x, p, n, 'ber')
% FORMAT ll  = spm_prob('Beta', 'logpdf', x, p, n, 'ber')
%   >> (Log) Probability density function.
%
% FORMAT e  = spm_prob('Beta', 'E',    p, n, 'ber')
% FORMAT el = spm_prob('Beta', 'Elog', p, n, 'ber')
% FORMAT v  = spm_prob('Beta', 'V',    p, n, 'ber')
% FORMAT vl = spm_prob('Beta', 'Vlog', p, n, 'ber')
%   >> Mean and variance.
%
% FORMAT kl  = spm_prob('Beta', 'kl', p1, n1, p0, n0, 'ber')
%   >> Kullback-Leibler divergence from B0 to B1 = KL(B1||B0).
%
% FORMAT [p1, n1] = spm_prob('Beta', 'up', p,  n,  p0, n0)
% FORMAT [p1, n1] = spm_prob('Beta', 'up', s0, s1, p0, n0, 'ber')
%   >> Posterior parameters of the Beta distribution.
%
%--------------------------------------------------------------------------
% Binomial probability conjugate
% ------------------------------
%
% To make the distribution easier to manipulate when used as a conjugate
% prior for a Binomial distribution, we reparameterise it by setting:
% - a = knp
% - b = kn(1-p)
% The expected value of the distribution becomes E[x] = p
% The variance of the distribution becomes       V[x] = p(1-p)/(kn+1)
% Then, the Beta distribution is parameterised by a mean probability  
% parameter (p), a degrees of freedom (n) and a number of trials (k).
%
% FORMAT pdf = spm_prob('Beta', 'pdf',    x, p, n, k, 'bin')
% FORMAT ll  = spm_prob('Beta', 'logpdf', x, p, n, k, 'bin')
%   >> (Log) Probability density function.
%
% FORMAT e  = spm_prob('Beta', 'E',    p, n, k, 'bin')
% FORMAT el = spm_prob('Beta', 'Elog', p, n, k, 'bin')
% FORMAT v  = spm_prob('Beta', 'V',    p, n, k, 'bin')
% FORMAT vl = spm_prob('Beta', 'Vlog', p, n, k, 'bin')
%   >> Mean and variance.
%
% FORMAT kl  = spm_prob('Beta', 'kl', p1, n1, p0, n0, k, 'bin')
%   >> Kullback-Leibler divergence from B0 to B1 = KL(B1||B0).
%
% FORMAT [p1, n1] = spm_prob('Beta', 'up', p,  n,  p0, n0)
% FORMAT [p1, n1] = spm_prob('Beta', 'up', s0, s1, p0, n0, k, 'bin')
%   >> Posterior parameters of the Beta distribution.
%
%--------------------------------------------------------------------------
% Negative-Binomial probability conjugate
% ---------------------------------------
%
% To make the distribution easier to manipulate when used as a conjugate
% prior for a negative Binomial distribution, we reparameterise it by 
% setting:
% - a = rnp/(1-p)
% - b = rn
% The expected value of the distribution becomes E[x] = p
% The variance of the distribution becomes       V[x] = p(1-p)^2/(rn+1-p)
% Then, the Beta distribution is parameterised by a mean probability  
% parameter (p), a degrees of freedom (n) and a number of failures (r).
%
% FORMAT pdf = spm_prob('Beta', 'pdf',    x, p, n, r, 'nbin')
% FORMAT ll  = spm_prob('Beta', 'logpdf', x, p, n, r, 'nbin')
%   >> (Log) Probability density function.
%
% FORMAT e  = spm_prob('Beta', 'E',    p, n, r, 'nbin')
% FORMAT el = spm_prob('Beta', 'Elog', p, n, r, 'nbin')
% FORMAT v  = spm_prob('Beta', 'V',    p, n, r, 'nbin')
% FORMAT vl = spm_prob('Beta', 'Vlog', p, n, r, 'nbin')
%   >> Mean and variance.
%
% FORMAT kl  = spm_prob('Beta', 'kl', p1, n1, p0, n0, r, 'nbin')
%   >> Kullback-Leibler divergence from B0 to B1 = KL(B1||B0).
%
% FORMAT [p1, n1] = spm_prob('Beta', 'up', p,  n,  p0, n0)
% FORMAT [p1, n1] = spm_prob('Beta', 'up', s0, s1, p0, n0, r, 'nbin')
%   >> Posterior parameters of the Beta distribution.
%
%--------------------------------------------------------------------------
% Geometric probability conjugate
% -------------------------------
%
% To make the distribution easier to manipulate when used as a conjugate
% prior for a Geometric distribution, we reparameterise it by 
% setting:
% - a = n
% - b = n(1-p)/p
% The expected value of the distribution becomes E[x] = p
% The variance of the distribution becomes       V[x] = (1-p)/(n+p)
% Then, the Beta distribution is parameterised by a mean probability  
% parameter (p) and a degrees of freedom (n).
%
% FORMAT pdf = spm_prob('Beta', 'pdf',    x, p, n, 'geom')
% FORMAT ll  = spm_prob('Beta', 'logpdf', x, p, n, 'geom')
%   >> (Log) Probability density function.
%
% FORMAT e  = spm_prob('Beta', 'E',    p, n, 'geom')
% FORMAT el = spm_prob('Beta', 'Elog', p, n, 'geom')
% FORMAT v  = spm_prob('Beta', 'V',    p, n, 'geom')
% FORMAT vl = spm_prob('Beta', 'Vlog', p, n, 'geom')
%   >> Mean and variance.
%
% FORMAT kl  = spm_prob('Beta', 'kl', p1, n1, p0, n0, 'geom')
%   >> Kullback-Leibler divergence from B0 to B1 = KL(B1||B0).
%
% FORMAT [p1, n1] = spm_prob('Beta', 'up', p,  n,  p0, n0)
% FORMAT [p1, n1] = spm_prob('Beta', 'up', s0, s1, p0, n0, 'geom')
% FORMAT [p1, n1] = spm_prob('Beta', 'up', s0, s1, p0, n0, 'fgeom')
%   >> Posterior parameters of the Beta distribution.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
    if nargin == 0
        help spm_prob>beta
        error('Not enough argument. Type ''help spm_prob>Beta'' for help.');
    end
    id = varargin{1};
    varargin = varargin(2:end);
    switch lower(id)
        case 'pdf'
            [varargout{1:nargout}] = beta_pdf(varargin{:});
        case {'logpdf', 'll'}
            [varargout{1:nargout}] = beta_logpdf(varargin{:});
        case 'kl'
            [varargout{1:nargout}] = beta_kl(varargin{:});
        case {'up', 'update'}
            [varargout{1:nargout}] = beta_up(varargin{:});
        case 'e'
            [varargout{1:nargout}] = beta_e(varargin{:});
        case 'elogdet'
            [varargout{1:nargout}] = beta_elog(varargin{:});
        case 'vloget'
            [varargout{1:nargout}] = beta_vlog(varargin{:});
        case 'help'
            help spm_prob>wishart
        otherwise
            help spm_prob>wishart
            error('Unknown function %s. Type ''help spm_prob>Wishart'' for help.', id)
    end
end

% -------------------------------------------------------------------------

function pdf = beta_pdf(x, varargin)
% FORMAT pdf = beta_pdf(x, a, b)
% FORMAT pdf = beta_pdf(x, p, n, 'ber')
% FORMAT pdf = beta_pdf(x, p, n, k, 'bin')
% FORMAT pdf = beta_pdf(x, p, n, r, 'nbin')
% FORMAT pdf = beta_pdf(x, p, n, 'geom')

    % Check if we are in the reparameterised case
    if nargin > 3
        switch lower(varargin{end})
            case 'ber'
                p = varargin{1};
                n = varargin{2};
                a = bsxfun(@times, n, p);
                b = bsxfun(@times, n, 1-p);
            case 'bin'
                p = varargin{1};
                n = varargin{2};
                k = varargin{3};
                a = bsxfun(@times, k, bsxfun(@times, n, p));
                b = bsxfun(@times, k, bsxfun(@times, n, 1-p));
            case 'nbin'
                p = varargin{1};
                n = varargin{2};
                r = varargin{3};
                b = bsxfun(@times, r, n);
                a = bsxfun(@rdivide, bsxfun(@times, p, b), 1-p);
            case 'geom'
                p = varargin{1};
                n = varargin{2};
                a = n;
                b = bsxfun(@rdivide, bsxfun(@times, n, 1-p), p);
        end
        pdf = beta_pdf(x, a, b);
        return
    end
    
    % Usual pdf
    a   = varargin{1};
    b   = varargin{2};
    pdf = betapdf(x, a, b);
end

% -------------------------------------------------------------------------

function pdf = beta_logpdf(x, varargin)
% FORMAT pdf = beta_logpdf(x, a, b)
% FORMAT pdf = beta_logpdf(x, p, n, 'ber')
% FORMAT pdf = beta_logpdf(x, p, n, k, 'bin')
% FORMAT pdf = beta_logpdf(x, p, n, r, 'nbin')
% FORMAT pdf = beta_logpdf(x, p, n, 'geom')

    % Check if we are in the reparameterised case
    if nargin > 3
        switch lower(varargin{end})
            case 'ber'
                p = varargin{1};
                n = varargin{2};
                a = bsxfun(@times, n, p);
                b = bsxfun(@times, n, 1-p);
            case 'bin'
                p = varargin{1};
                n = varargin{2};
                k = varargin{3};
                a = bsxfun(@times, k, bsxfun(@times, n, p));
                b = bsxfun(@times, k, bsxfun(@times, n, 1-p));
            case 'nbin'
                p = varargin{1};
                n = varargin{2};
                r = varargin{3};
                b = bsxfun(@times, r, n);
                a = bsxfun(@rdivide, bsxfun(@times, p, b), 1-p);
            case 'geom'
                p = varargin{1};
                n = varargin{2};
                a = n;
                b = bsxfun(@rdivide, bsxfun(@times, n, 1-p), p);
        end
        pdf = beta_logpdf(x, a, b);
        return
    end
    
    % Usual log-pdf
    a   = varargin{1};
    b   = varargin{2};
    pdf = bsxfun(@plus,  bsxfun(@times, (a-1), log(x)), ...
                         bsxfun(@times, (b-1), log(1-x)));
    pdf = bsxfun(@minus, pdf, betaln(a, b));
end

% -------------------------------------------------------------------------

function e = beta_e(varargin)
% FORMAT e = beta_e(a, b)
% FORMAT e = beta_e(p, n, 'ber')
% FORMAT e = beta_e(p, n, k, 'bin')
% FORMAT e = beta_e(p, n, r, 'nbin')
% FORMAT e = beta_e(p, n, 'geom')

    % Check if we are in the reparameterised case
    if nargin > 2
        p = varargin{1};
        e = p;
        return
    end
    
    % Usual expected value
    % E[x] = a/(a+b)
    a = varargin{1};
    b = varargin{2};
    e = bsxfun(@rdivide, a, bsxfun(@plus, a, b));
end

% -------------------------------------------------------------------------

function el = beta_elog(varargin)
% FORMAT el = beta_elog(a, b)
% FORMAT el = beta_elog(p, n, 'ber')
% FORMAT el = beta_elog(p, n, k, 'bin')
% FORMAT el = beta_elog(p, n, r, 'nbin')
% FORMAT el = beta_elog(p, n, 'geom')

    % Check if we are in the reparameterised case
    if nargin > 2
        switch lower(varargin{end})
            case 'ber'
                p = varargin{1};
                n = varargin{2};
                a = bsxfun(@times, n, p);
                b = bsxfun(@times, n, 1-p);
            case 'bin'
                p = varargin{1};
                n = varargin{2};
                k = varargin{3};
                a = bsxfun(@times, k, bsxfun(@times, n, p));
                b = bsxfun(@times, k, bsxfun(@times, n, 1-p));
            case 'nbin'
                p = varargin{1};
                n = varargin{2};
                r = varargin{3};
                b = bsxfun(@times, r, n);
                a = bsxfun(@rdivide, bsxfun(@times, p, b), 1-p);
            case 'geom'
                p = varargin{1};
                n = varargin{2};
                a = n;
                b = bsxfun(@rdivide, bsxfun(@times, n, 1-p), p);
        end
        el = beta_elog(a, b);
        return
    end
    
    % Usual expected value
    % E[ln x] = psi(a) - psi(a+b)
    a  = varargin{1};
    b  = varargin{2};
    el = bsxfun(@minus, psi(a), psi(bsxfun(@plus, a, b)));
end

% -------------------------------------------------------------------------

function v = beta_v(varargin)
% FORMAT v = beta_v(a, b)
% FORMAT v = beta_v(p, n, 'ber')
% FORMAT v = beta_v(p, n, k, 'bin')
% FORMAT v = beta_v(p, n, r, 'nbin')
% FORMAT v = beta_v(p, n, 'geom')

    % Check if we are in the reparameterised case
    if nargin > 2
        switch lower(varargin{end})
            case 'ber'
                p = varargin{1};
                n = varargin{2};
                a = bsxfun(@times, n, p);
                b = bsxfun(@times, n, 1-p);
            case 'bin'
                p = varargin{1};
                n = varargin{2};
                k = varargin{3};
                a = bsxfun(@times, k, bsxfun(@times, n, p));
                b = bsxfun(@times, k, bsxfun(@times, n, 1-p));
            case 'nbin'
                p = varargin{1};
                n = varargin{2};
                r = varargin{3};
                b = bsxfun(@times, r, n);
                a = bsxfun(@rdivide, bsxfun(@times, p, b), 1-p);
            case 'geom'
                p = varargin{1};
                n = varargin{2};
                a = n;
                b = bsxfun(@rdivide, bsxfun(@times, n, 1-p), p);
        end
        v = beta_v(a, b);
        return
    end
    
    % Usual variance
    % V[x] = ab/[(a+b)^2*(a+b+1)]
    a  = varargin{1};
    b  = varargin{2};
    v = bsxfun(@plus, a, b);
    v = bsxfun(@times, v.^2, v+1);
    v = bsxfun(@rdivide, bsxfun(@times, a, b), v);
end

% -------------------------------------------------------------------------

function vl = beta_vlog(varargin)
% FORMAT vl = beta_vlog(a, b)
% FORMAT vl = beta_vlog(p, n, 'ber')
% FORMAT vl = beta_vlog(p, n, k, 'bin')
% FORMAT vl = beta_vlog(p, n, r, 'nbin')
% FORMAT vl = beta_vlog(p, n, 'geom')

    % Check if we are in the reparameterised case
    if nargin > 2
        switch lower(varargin{end})
            case 'ber'
                p = varargin{1};
                n = varargin{2};
                a = bsxfun(@times, n, p);
                b = bsxfun(@times, n, 1-p);
            case 'bin'
                p = varargin{1};
                n = varargin{2};
                k = varargin{3};
                a = bsxfun(@times, k, bsxfun(@times, n, p));
                b = bsxfun(@times, k, bsxfun(@times, n, 1-p));
            case 'nbin'
                p = varargin{1};
                n = varargin{2};
                r = varargin{3};
                b = bsxfun(@times, r, n);
                a = bsxfun(@rdivide, bsxfun(@times, p, b), 1-p);
            case 'geom'
                p = varargin{1};
                n = varargin{2};
                a = n;
                b = bsxfun(@rdivide, bsxfun(@times, n, 1-p), p);
        end
        vl = beta_vlog(a, b);
        return
    end
    
    % Usual variance
    % V[ln x] = psi_1(a) - psi_1(a+b)
    a  = varargin{1};
    b  = varargin{2};
    vl = bsxfun(@minus, psi(1,a), psi(1,bsxfun(@plus, a, b)));
end

% -------------------------------------------------------------------------

function kl = beta_kl(varargin)
% FORMAT vl = beta_kl(a1, b1, a0, b0)
% FORMAT vl = beta_kl(p1, n1, p0, n0, 'ber')
% FORMAT vl = beta_kl(p1, n1, p0, n0, k, 'bin')
% FORMAT vl = beta_kl(p1, n1, p0, n0, r, 'nbin')
% FORMAT vl = beta_kl(p1, n1, 'geom')

    % Check if we are in the reparameterised case
    if nargin > 4
        switch lower(varargin{end})
            case 'ber'
                p1 = varargin{1};
                n1 = varargin{2};
                p0 = varargin{3};
                n0 = varargin{4};
                a1 = bsxfun(@times, n1, p1);
                b1 = bsxfun(@times, n1, 1-p1);
                a0 = bsxfun(@times, n0, p0);
                b0 = bsxfun(@times, n0, 1-p0);
            case 'bin'
                p1 = varargin{1};
                n1 = varargin{2};
                p0 = varargin{3};
                n0 = varargin{4};
                k  = varargin{5};
                a1 = bsxfun(@times, k, bsxfun(@times, n1, p1));
                b1 = bsxfun(@times, k, bsxfun(@times, n1, 1-p1));
                a0 = bsxfun(@times, k, bsxfun(@times, n0, p0));
                b0 = bsxfun(@times, k, bsxfun(@times, n0, 1-p0));
            case 'nbin'
                p1 = varargin{1};
                n1 = varargin{2};
                p0 = varargin{3};
                n0 = varargin{4};
                r  = varargin{5};
                b1 = bsxfun(@times, r, n1);
                a1 = bsxfun(@rdivide, bsxfun(@times, p1, b1), 1-p1);
                b0 = bsxfun(@times, r, n0);
                a0 = bsxfun(@rdivide, bsxfun(@times, p0, b0), 1-p0);
            case 'geom'
                p1 = varargin{1};
                n1 = varargin{2};
                p0 = varargin{3};
                n0 = varargin{4};
                a1 = n1;
                b1 = bsxfun(@rdivide, bsxfun(@times, n1, 1-p1), p1);
                a0 = n0;
                b0 = bsxfun(@rdivide, bsxfun(@times, n0, 1-p0), p0);
        end
        kl = beta_kl(a1, b1, a0, b0);
        return
    end
    
    % Usual KL-divergence
    % KL(a1,b1||a0,b0) = ln B(a0,b0) - ln B(a1,b1) 
    %                    + (a1-a0)psi(a1) + (b1-b0)psi(b1)
    %                    + (a0-a1+b0-b1)psi(a1+b1)
    a1 = varargin{1};
    b1 = varargin{2};
    a0 = varargin{3};
    b0 = varargin{4};
    kl = bsxfun(@times, ...
        bsxfun(@plus, bsxfun(@minus, a0, a1), bsxfun(@minus, b0, b1)), ...
        psi(bsxfun(@plus, a1, b1)));
    kl = bsxfun(@plus,  kl, bsxfun(@times, bsxfun(@minus, a1, a0), psi(a1)));
    kl = bsxfun(@plus,  kl, bsxfun(@times, bsxfun(@minus, b1, b0), psi(b1)));
    kl = bsxfun(@plus,  kl, betaln(a0,b0));
    kl = bsxfun(@minus, kl, betaln(a1,b1));
end

% -------------------------------------------------------------------------

function [p1,n1] = beta_up(varargin)
% FORMAT [p1,n1] = beta_up(p, n, p0, n0)
%
% FORMAT [p1,n1] = beta_up(s0, s1, p0, n0,    'suffstat', 'ber')   p=s1/s0
% FORMAT [p1,n1] = beta_up(s0, s1, p0, n0, k, 'suffstat', 'bin')   p=s1/(k*s0)
% FORMAT [p1,n1] = beta_up(s0, s1, p0, n0, r, 'suffstat', 'nbin')  p=s1/(r*s0+s1)
% FORMAT [p1,n1] = beta_up(s0, s1, p0, n0,    'suffstat', 'geom')  p=s0/s1
% FORMAT [p1,n1] = beta_up(s0, s1, p0, n0,    'suffstat', 'fgeom') p=s0/(s1+s0)

    if ischar(varargin{end})
        switch lower(varargin{end})
            case 'ber'
                s0 = varargin{1};
                s1 = varargin{2};
                n = s0;
                p = s1/s0;
            case 'bin'
                s0 = varargin{1};
                s1 = varargin{2};
                k  = varargin{5};
                n = s0;
                p = s1/(k*s0);
            case 'nbin'
                s0 = varargin{1};
                s1 = varargin{2};
                r  = varargin{5};
                n = s0;
                p = s1/(r*s0+s1);
            case 'geom'
                s0 = varargin{1};
                s1 = varargin{2};
                n = s0;
                p = s0/s1;
            case 'fgeom'
                s0 = varargin{1};
                s1 = varargin{2};
                n = s0;
                p = s0/(s1+s0);
        end
    else
        p = varargin{1};
        n = varargin{2};
    end
    p0 = varargin{3};
    n0 = varargin{4};
    n1 = n0+n;
    p1 = (n0*p0+n*p)/n1;

end
