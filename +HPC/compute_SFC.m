% main SFC function

% modular version to run on HPC
function compute_SFC(cpu_nr, n_iter, cellArea, lfpArea, cds, chanType, cellType)

Recall.HPC.setPaths_sfc_hpc
rng('shuffle')
if nargin == 6, cellType = 'sigRamp'; end % might change later

% % should just always use this 
% if strcmp(cellArea(2:end), 'FFA') && isempty(cellType)
%    cellType = 'sigRamp';
% end

if strcmp(cds, 'ScreeningImagination')
    shortcds = 'SI';
elseif strcmp(cds, 'EncodingImagination')
    shortcds = 'EI';
end

fname = [diskPath filesep 'Recall_Task' filesep 'ppc_log' filesep 'Data' filesep ['data_SFC_' cellArea 'Cell_' lfpArea 'LFP' '_' cellType '_' cds '_' chanType]];

% if the matfile actually exists
assert(exist([fname '.mat']) == 2, 'File not found - check name');
load(fname)

n_freq = params.high_freq - params.low_freq;
for sess = 1:size(data_lfp, 1)
    if ~isempty(data_lfp{sess, 1}) && ~isempty(data_lfp{sess, 2})
        if strcmp(params.balance_spikes, 'true')
            
            for n = 1:n_iter % n_desiredIter/n_cores assigned
                
                [d1, d2, resampledCond] = Utilities.LFP.balanceSpikesSFC(data_spike{sess, 1}, data_spike{sess, 2});
                assert(isequal(sum(d1(:)), sum(d2(:))));
                
                
                % compute ppc for that session and condition
                if resampledCond == 1
                    downSampledData_spike = d1;
                    unchangedData_spike = d2;
                    unchangedCond = 2;
                elseif resampledCond == 2
                    downSampledData_spike = d2;
                    unchangedData_spike = d1;
                    unchangedCond = 1;
                end
                
                if ~isempty(data_lfp{sess, resampledCond})
                    [ppc{sess, resampledCond}(:, :, :, n), frq, ppc_boot{sess, resampledCond}(:, :, :, :, n)]...
                        = Recall.compute_ppc(data_lfp{sess, resampledCond}, downSampledData_spike, params.low_freq, params.high_freq, n_freq, params.scale, params.FsDown, params.run_boot);
                end
            end
            
            % run the unchanged condition just once
            if ~isempty(data_lfp{sess, unchangedCond})
                [ppc{sess, unchangedCond}, frq, ppc_boot{sess, unchangedCond}]...
                    = Recall.compute_ppc(data_lfp{sess, unchangedCond}, unchangedData_spike, params.low_freq, params.high_freq, n_freq, params.scale, params.FsDown, params.run_boot);
            end
        else
            n_iter = 0;
            for cond = 1:2
                % run the unchanged condition just once
                if ~isempty(data_lfp{sess, unchangedCond})
                    [ppc{sess, cond}, frq, ppc_boot{sess, cond}]...
                        = Recall.compute_ppc(data_lfp{sess, cond}, data_spike{sess, cond}, params.low_freq, params.high_freq, n_freq, params.scale, params.FsDown, params.run_boot);
                end
            end
        end
    end
end

% Save 
outPath = [diskPath filesep 'Recall_Task' filesep 'ppc_' params.scale];
outDir = [outPath filesep ['ppc_' params.cellArea 'Cell_' params.lfpArea 'LFP_' params.cells '_' shortcds  '_' chanType]];
if ~exist(outDir, 'dir')
    mkdir(outDir)
end

fnum = sprintf('%03d', cpu_nr);
filename = [outDir filesep ['ppc_' params.cellArea '_' params.lfpArea '_' shortcds '_' chanType '_subsampleIterations_' num2str(n_iter) '_worker_' fnum]];
% disp(filename)
save(filename, 'ppc', 'ppc_boot', 'frq', 'params', 'fname');



end