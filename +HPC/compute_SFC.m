% main SFC script

% modular versio nto run on HPC

setPaths_sfc_hpc


%% Execute

tic   
for sess = 1:length(dirID)
    
    if sess ~= 1
        params = rmfield(params, 'valid_cell');
    end
    
    sessDir = dirID(sess, :);
    
    [lfpChans] = Utilities.LFP.defineLFPChannelsSFC(sessDir, lfpArea);
    
    % read in session IDs and assign conditions etc.
    [params] = Utilities.LFP.defineInputParamsSFC(cellArea, lfpArea, lfpChans, sessDir{1}, 'log', 'sigRamp', diskPath);
    
    for cond = 1:2
        
        % note the directories are different for screening vs Im and Encoding vs Im
        if cond == 1
            condition = condition_1;
            
        elseif cond == 2
            condition = condition_2;
            
        end
     
        % these structs contain time windows for lfp and spikes
        % other logistical things needed for data extraction
        [lfDat, spikDat] = Utilities.LFP.SFCConfig(params, condition);
        
        % grab data - including choosing valid data
        [data_lfp{sess, cond}, data_spike{sess, cond}, params] = Utilities.LFP.ExtractDataSFC(lfDat, spikDat, params, sessDir, condition);
        
        % compute SFC for valid cell/channel pairs
        n_freq = params.high_freq - params.low_freq;
   
    end
end
toc
 

%%
for sess = 1:length(dirID)-1
        
        
    if strcmp(params.balance_spikes, 'true')
        for n = 1:200
            
            [d1, d2] = Utilities.LFP.balanceSpikesSFC(data_spike{sess, 1}, data_spike{sess, 2});
            
            assert(isequal(sum(d1(:)), sum(d2(:))));
            
            %             % compute ppc for that session and condition
            %             if ~isempty(data_lfp)
            %                 [ppc{sess, cond}, frq, ppc_boot{sess, cond}] = Recall.compute_ppc(data_lfp{sess, cond}, data_spike{sess, cond}, params.low_freq, params.high_freq, n_freq, params.scale, params.FsDown, params.run_boot);
            %             end
        end
    end
    
end

%% Save 
if strcmp(params.scale, 'log')
    save([diskPath filesep 'Recall_task' filesep 'ppc_log' filesep ['ppc_' params.cellArea 'Cell_' params.lfpArea 'LFP_' params.cells '_' cds  '_cellChans']], 'ppc', 'ppc_boot', 'frq', 'params');
elseif strcmp(params.scale, 'linear')
    save([diskPath filesep 'Recall_task' filesep 'ppc_linear' filesep ['ppc_' params.cellArea 'Cell_' params.lfpArea 'LFP_' params.cells '_' cds  '_cellChans']], 'ppc', 'ppc_boot', 'frq', 'params');
end

%% Go to compute_SFC_Stats or plot_SFC_main 



