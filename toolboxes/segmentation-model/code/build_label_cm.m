function [dat,CM,opt] = build_label_cm(dat,opt)
% FORMAT [dat,CM,opt] = build_label_cm(dat,opt)
% dat - Subjects data structure
% opt - Options structure
% CM  - Confusion matrix
%
% Build Rater confusion matrix for each subject.
% This matrix maps template classes to manually segmented classes.
% Manual labels often do not follow the same convention as the Template, 
% and not all regions may be labelled. Therefore, a manual label may 
% correspond to several Template classes and, conversely, one Template
% class may correspond to several manual labels.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Parameters
S0 = numel(dat);

% For each subject with labels, build the label confusion matrix
for s=1:S0        
    
    CM = get_label_cm(dat{s}.population,opt);
    
    % Assign to dat
    dat{s}.gmm.cm = CM;
end
%===========================================================================