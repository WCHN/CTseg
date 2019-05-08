function out = ctseg_run(job)

if isempty(job.outdir{1})  
    outdir = '.';   
else
    outdir = job.outdir{1}; 
end
if isempty(job.cleanbrain) 
    cleanbrain = false; 
else
    cleanbrain = job.cleanbrain; 
end
if isempty(job.image) 
    image = [1 0]; 
else
    image = job.image; 
end
if isempty(job.native) 
    native = [1 0]; 
else
    native = job.native; 
end
if isempty(job.warped) 
    warped = [0 0]; 
else
    warped = job.warped; 
end
if isempty(job.samp) 
    Samp = 2; 
else
    Samp = job.samp; 
end
if isempty(job.mrf) 
    MRF = 1; 
else
    MRF = job.mrf; 
end

Write        = struct;
Write.image  = image;
Write.native = native;
Write.warped = warped;

N       = size(job.data,1);
results = cell(1,N);
for n=1:N
    Image      = deblank(job.data{n});    
    results{n} = spm_segment_ct(Image,outdir,cleanbrain,Write,Samp,MRF);
end

% TODO: make possible to use dependecies
out = [];
% tiss = struct('c',{});
% for j=1:8
%     tiss(j).c = cell(N,1);
% end
% for n=1:N
%     [pth,nam,ext] = fileparts(job.data{n});
%     for j=1:8
%         tiss(j).c{n} = fullfile(outdir,['c',num2str(j),nam,'.nii']);
%     end
% end
% out.tiss = tiss;
