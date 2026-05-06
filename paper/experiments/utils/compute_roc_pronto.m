function [fpr, tpr, auc] = compute_roc_pronto(PRT, model)
% Compute ROC curve from a PRoNTo PRT struct, using func_val (continuous
% GP posterior probabilities). Follows prt_plot_ROC.m logic exactly.
%
% FORMAT [fpr, tpr, auc] = compute_roc_pronto(PRT, model)
%
% INPUT
%   PRT   - PRoNTo data structure with at least one estimated model
%   model - model number (usually 1)
%
% OUTPUT
%   fpr - false positive rates
%   tpr - true positive rates
%   auc - area under the ROC curve

    % Collect func_val and targets across all folds
    nfold   = length(PRT.model(model).output.fold);
    fVals   = [];
    targets = [];
    for f = 1:nfold
        targets = [targets; PRT.model(model).output.fold(f).targets]; %#ok<AGROW>
        if isfield(PRT.model(model).output.fold(f), 'func_val')
            fVals = [fVals; PRT.model(model).output.fold(f).func_val]; %#ok<AGROW>
        else
            fVals = [fVals; PRT.model(model).output.fold(f).predictions]; %#ok<AGROW>
        end
    end

    % Compute TPR and FPR by sweeping sorted thresholds (as in prt_plot_ROC)
    targpos = targets == 1;
    numPos  = sum(targpos);
    numNeg  = sum(~targpos);

    [s_scores, idx] = sort(fVals, 'descend');
    s_targets = targets(idx);

    tp = 0; fp = 0;
    thr_prev   = s_scores(1) + 1;
    num_scores = length(s_scores);
    num_thr    = length(unique(s_scores)) + 1;
    tpr = []; fpr = [];
    i = 1; j = 1;
    while i <= num_scores && j <= num_thr
        if thr_prev ~= s_scores(i)
            tpr(j,1) = tp / numPos;
            fpr(j,1) = fp / numNeg;
            thr_prev = s_scores(i);
            j = j + 1;
        end
        if s_targets(i) == 1
            tp = tp + 1;
        else
            fp = fp + 1;
        end
        i = i + 1;
    end
    tpr(j,1) = tp / numPos;
    fpr(j,1) = fp / numNeg;

    % AUC via trapezoidal rule
    n   = size(tpr, 1);
    auc = sum((fpr(2:n) - fpr(1:n-1)) .* (tpr(2:n) + tpr(1:n-1))) / 2;
end
