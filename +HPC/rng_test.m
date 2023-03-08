function rng_test(cpu_nr, testInput)
% Quick function to test random number generation
% on cluster and a first slurm script
% vwadia March2023
Recall.HPC.setPaths_sfc_hpc;

r = randi(100, [5 1]);

outDir = [diskPath filesep 'Recall_Task' filesep 'ppc_log' filesep 'Data'];
if ~exist(outDir, 'dir')
    mkdir(outDir)
end
save([outDir filesep ['rng_worker_' num2str(cpu_nr)]], 'r', 'testInput');