% version to run on HPC
% computing ssim on all things takes forever (~6 months)
function compute_image_similarity(cpu_nr, n_image_start, n_image_end)

Screening.HPC.setPaths_imageSim_hpc
rng('shuffle')


% easier to just have all paths point to this folder
fname = [diskPath filesep 'Recall_Task' filesep 'Imagecd ../Data' filesep 'LargeImageSet_17856ims'];
disp(fname);

% if the matfile actually exists
assert(exist([fname '.mat']) == 2, 'File not found - check name');
load(fname)

imageIDs = [1:500]';
bigSetIds = 1:size(im2all_woStim, 3);

ssimVals_imIDs = cell(length(imageIDs), 1);


for ctr = n_image_start:n_image_end
    targ = im2all_woStim(:, :, ctr);
%     tic
    for ctr2 = ctr:size(im2all_woStim, 3)
        comp = im2all_woStim(:, :, ctr2);
        simVal(ctr2) = ssim(targ, comp); 
        
    end
    ssimVals_imIDs{ctr} = simVal;
%     toc
end

% Save
outPath = [diskPath filesep 'Recall_Task' filesep 'ImageData'];
if ~exist(outPath, 'dir')
    mkdir(outPath)
end


fnum = sprintf('%03d', cpu_nr);
filename = [outPath filesep 'SSIMVals_Images_' num2str(n_image) '_' fnum];
disp(filename)
save(filename, 'ssimVals_imIDs');

end