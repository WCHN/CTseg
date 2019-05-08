function modality = get_modality_name(dat,name_population)
S = numel(dat);
for s=1:S
    population = dat{s}.population;

    if strcmp(name_population,population)
        
        if isfield(dat{s}.modality{1},'channel') && numel(dat{s}.modality{1}.channel) > 1
            modality = '';
        elseif isfield(dat{s}.modality{1},'channel')
            modality = dat{s}.modality{1}.channel{1}.name;
        else
            modality = dat{s}.modality{1}.name;
        end
        break
    end
end
%==========================================================================