function [icc, lb, ub] = compute_icc(x, y, alpha)
% Intraclass correlation coefficient ICC(3,1) — two-way mixed, single measures.
% Used for assessing agreement between two methods measuring the same quantity.
%
%   x, y:   column vectors of measurements from two methods
%   alpha:  significance level for CI (default: 0.05)
%   icc:    ICC(3,1) value
%   lb, ub: lower and upper bounds of 95% CI

    if nargin < 3, alpha = 0.05; end

    n = numel(x);
    data = [x(:), y(:)];
    k = 2;  % number of raters/methods

    % Grand mean
    gm = mean(data(:));

    % Sum of squares
    row_means = mean(data, 2);
    col_means = mean(data, 1);

    SSR = k * sum((row_means - gm).^2);          % between subjects
    SSC = n * sum((col_means - gm).^2);          % between methods
    SST = sum((data(:) - gm).^2);                % total
    SSE = SST - SSR - SSC;                        % residual

    % Mean squares
    MSR = SSR / (n - 1);
    MSE = SSE / ((n - 1) * (k - 1));
    MSC = SSC / (k - 1);

    % ICC(3,1)
    icc = (MSR - MSE) / (MSR + (k - 1) * MSE);

    % F-test for CI
    F = MSR / MSE;
    df1 = n - 1;
    df2 = (n - 1) * (k - 1);

    F_lo = F / finv(1 - alpha/2, df1, df2);
    F_hi = F * finv(1 - alpha/2, df2, df1);

    lb = (F_lo - 1) / (F_lo + k - 1);
    ub = (F_hi - 1) / (F_hi + k - 1);
end
