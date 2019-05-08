function dat = init_bf(dat,opt,scl)
% FORMAT dat = init_bf(dat,opt,[scl])
% dat   - Subjects data structure
% opt   - Options structure
% scl   - Scaling of the data (that kind of aligns histograms)
%
% Init registration related variables:
% * dat.bf.chan(:).C:       Bias precision matrix
% * dat.bf.chan(:).B1/2/3:  DCT basis functions along each dimension
% * dat.bf.chan(:).T:       Bias coordinates in DCT basis (all zero)
%
% The histogram scaling is converted to a component of the bias field.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if opt.bf.do
    biasreg = opt.bf.biasreg;
    fwhm    = opt.bf.biasfwhm;
    S0      = numel(dat);

    % Loop over all subjects
    for s=1:S0    
        modality = dat{s}.modality{1}.name;
        
%         if strcmpi(modality,'MRI')
            [dm,~,vs,C] = obs_info(dat{s});
            ff          = get_ff(vs);               

            if nargin < 3
                [~,~,~,~,scl1] = get_obs(dat{s},'mskonlynan',opt.seg.mskonlynan); % Do not subsample!
                scl{s}         = scl1;
            end

            [~,grd] = get_subsampling_grid(dm,vs,opt.seg.samp);

            cl             = cell(C,1);    
            args           = {'C',cl,'B1',cl,'B2',cl,'B3',cl,'T',cl,'ll',cl};
            dat{s}.bf.chan = struct(args{:});

            sd = vs(1)*dm(1)/fwhm; d3(1) = ceil(sd*2);
            sd = vs(2)*dm(2)/fwhm; d3(2) = ceil(sd*2);
            sd = vs(3)*dm(3)/fwhm; d3(3) = ceil(sd*2);

            prec  = spm_bias_lib('regulariser','bending',dm,d3,vs);
            prec  = prec*biasreg*ff;

            for c=1:C
                % GAUSSIAN REGULARISATION for bias correction                        
                dat{s}.bf.chan(c).C = prec;

                % Basis functions for bias correction
                dat{s}.bf.chan(c).B3  = spm_dctmtx(dm(3),d3(3),grd.z0);
                dat{s}.bf.chan(c).B2  = spm_dctmtx(dm(2),d3(2),grd.y0(1,:)');
                dat{s}.bf.chan(c).B1  = spm_dctmtx(dm(1),d3(1),grd.x0(:,1));

                % Initial parameterisation of bias field
                dat{s}.bf.chan(c).T = zeros(d3);

                % Change DC component of bias field to make intensities more
                % simillar between MR images.
                b1 = dat{s}.bf.chan(c).B1(1,1);
                b2 = dat{s}.bf.chan(c).B2(1,1);
                b3 = dat{s}.bf.chan(c).B3(1,1);

                dat{s}.bf.chan(c).T(1,1,1) = 1/(b1*b2*b3)*log(scl{s}(c));
            end
%         end
    end
end
%==========================================================================