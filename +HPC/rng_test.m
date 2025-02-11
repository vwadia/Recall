function rng_test(cpu_nr, testInput)
% Quick function to test random number generation
% on cluster and a first slurm script
% vwadia March2023
Recall.HPC.setPaths_sfc_hpc;
rng('shuffle')

r = randi(100, [5 1]);

outDir = [diskPath filesep 'Recall_Task' filesep 'ppc_log' filesep 'Data'];
if ~exist(outDir, 'dir')
    mkdir(outDir)
end
fnum = sprintf('%03d', cpu_nr);
save([outDir filesep ['rng_worker_' fnum]], 'r', 'testInput');