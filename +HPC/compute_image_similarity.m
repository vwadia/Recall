% version to run on HPC
% computing ssim on all things takes forever (~6 months)
function compute_image_similarity(cpu_nr)

Screening.HPC.setPaths_imageSim_hpc
rng('shuffle')


% easier to just have all paths point to this folder
fname = [diskPath filesep 'Recall_Task' filesep 'ImageData' filesep 'LargeImageSet_17856ims'];
disp(fname);

% if the matfile actually exists
assert(exist([fname '.mat']) == 2, 'File not found - check name');
load(fname)

imageIDs = [1:500]';
bigSetIds = 1:size(im2all_woStim, 3);

ssimVals_imIDs = cell(length(imageIDs), 1);

% hard coded numbers - im set is 17856, # of workers is 31 and 17856/31 = 576
n_image_start = (576*(cpu_nr-1))+1; 
n_image_end = 576*cpu_nr;

for ctr = n_image_start:n_image_end
    targ = im2all_woStim(:, :, ctr);
%     tic
    simVal = nan(1, size(im2all_woStim, 3));
    for ctr2 = 1:size(im2all_woStim, 3)
        comp = im2all_woStim(:, :, ctr2);
        simVal(ctr2) = ssim(targ, comp); 
        
    end
    ssimVals_imIDs{ctr} = simVal;
%     toc
end

% Save
outPath = [diskPath filesep 'Recall_Task' filesep 'ImageData' filesep 'hpc_output'];
if ~exist(outPath, 'dir')
    mkdir(outPath)
end


fnum = sprintf('%03d', cpu_nr);
filename = [outPath filesep 'SSIMVals_Images_' num2str(n_image_start) '_to_' num2str(n_image_start) '_worker_' fnum];
disp(filename)
save(filename, 'ssimVals_imIDs');

end