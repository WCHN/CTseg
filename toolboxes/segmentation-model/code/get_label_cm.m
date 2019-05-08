function CM = get_label_cm(dat,opt)
% FORMAT CM = get_label_cm(dat,opt)
% dat - Subjects data structure
% opt - Options structure
% CM  - confusion matrix
%
% Build Rater confusion matrix for one subject.
% This matrix maps template classes to manually segmented classes.
% Manual labels often do not follow the same convention as the Template, 
% and not all regions may be labelled. Therefore, a manual label may 
% correspond to several Template classes and, conversely, one Template
% class may correspond to several manual labels.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Here, we assume that all subjects from the same population (e.g.,
% a publicily available dataset) have the same labelling protocole and 
% confusion matrix.
% We allow the rater's sensitivity to change every few acquistion. We would
% typically start with a high sensitivity, to weight the labels strongly,
% and then decrease this value to allow the model to correct the rater's
% mistakes (especially near boundaries).

population = dat.population;    % Population of the subject
iter       = opt.model.it;      % Current iteration
cm         = opt.gmm.labels.cm; % Dictionary pop -> class mapping
rater_sens = opt.sched.labels(min(numel(opt.sched.labels),iter));
K          = opt.template.K;    % Number of template classes

if ~cm.isKey(population)
    AL = zeros(1,K);
else
    AL = cm(population);
end

%--------------------------------------------------------------------------
% Build confusion matrix 
%--------------------------------------------------------------------------

%% 
% M = [0 0 0 0 1 0 0 0 0 0; ...
%      0 0 0 0 0 0 0 1 0 0; ...
%      0 0 0 0 0 0 1 0 0 0; ...
%      0 0 1 0 0 0 0 0 0 0; ...
%      0 0 0 0 0 0 1 1 0 0; ...
%      0 0 0 0 0 0 0 0 1 0; ...
%      1 1 0 1 0 1 0 0 0 1];
% 
% M = {[5],[8],[7],[3],[7 8],[9],{[1 2 4 5 10],[9]}};

L  = numel(AL);    % Number of labels
CM = zeros([L K]); % Allocate confusion matrix
for l=1:L % Loop over labels
    
    ix        = false(1,K);
    ix(AL{l}) = true;
    
    CM(l,ix)  = rater_sens/nnz(ix); 
    CM(l,~ix) = (1 - rater_sens)/nnz(~ix);
end

CM = bsxfun(@rdivide,CM,sum(CM,2));

% ix_ul = max(ix) + 1;    % Unlabelled voxels (0) are assumed to have value: num(labels) + 1
% CM    = zeros(ix_ul,K); % Confusion matrix is of size: (num(labels) + 1) x num(tissue)
% ix_k  = ix == 0; 
% 
% % For the unlabelled voxels, we assume close to uniform probability
% CM(ix_ul, ix_k) = opt.gmm.labels.Su;
% CM(ix_ul,~ix_k) = 1 - opt.gmm.labels.Su;    
% 
% % For the labelled voxels, we define probabilities from a user-given
% % rater sensitivity
% for k=ix
%     if k == 0
%         % Skip unlabelled
%         continue; 
%     end
% 
%     ix_k        = ix == k;
%     CM(k, ix_k) = rater_sens/nnz(ix_k);
%     CM(k,~ix_k) = (1 - rater_sens)/nnz(~ix_k);
% end
% 
% % Normalise confusion matrix
% CM = bsxfun(@rdivide,CM,sum(CM,2));

return