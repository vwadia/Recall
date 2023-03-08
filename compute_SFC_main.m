% main SFC script

% Flow:
%     1. DefineInputParams (run as script not function?)
%         a. direction (A cell to B lfp and vice versa)
%         b. condition (task condition)
%         c. session
%         note: Just do L and R in all cases
%         Maybe just do both directions in all cases (or if not defined)?
%         
%     2. SFCConfig - define session and read in input parameters 
%     
%     3. Grab the spike rasters and the lfp snippets for the relevant condition, brain areas etc. 
%     
%     4. DefineValidData - big endeavour
%         For Spikes:
%             a. Canonical way: count spikes in baseline period (require more than 50)  
%             b. Could also pre-select visually responsive cells        
%         For Channels:
%             a. Canonical way: only keep channels that had spikes on them
%         Cell-Channel Pairs:
%             a. Could compute surrogate distributions:
%                 Take spike times and lfp trial snippets - compute sfc
%                 Then take spike times and SHUFFLED lfp trial snippets - compute sfc again
%                 is the former significantly different from the latter? If yes then choose channel
%             b. Best to compute these distributions using trials from all conditions
%                 Or equally sample channels that have sig sfc in either condition
%                 
%     5. Compute ppc for all cell-channel pairs remaining

% vwadia/Jan2023

%% set paths and define parameters

setDiskPaths
taskCodePath = [boxPath filesep 'recallTaskVarun']; % more events described here
addpath(taskCodePath); setTTLCodes;
addpath([diskPath filesep 'Code' filesep 'SFCpackage' filesep 'helpers']);
addpath([diskPath filesep 'Code' filesep 'SFCpackage' filesep 'SFCieeg']);
addpath(genpath([diskPath filesep 'Code' filesep 'osortTextUI']));

saveDataOnly = 1; % for cluster runs

%% Define patients, areas, and channel list 

patientIDs = {'P76CS', 'P79CS', 'P80CS'};
cellArea = 'LFFA'; lfpArea = 'LH';
% cellArea = 'LH'; lfpArea = 'LFFA';

condition_1 = 'Screening'; cds = ['ScreeningImagination'];
% condition_1 = 'Encoding'; cds = ['EncodingImagination'];
condition_2 = 'Imagination';

chanType = 'Cell';;
% chanType = 'NoNoise';
dirID = Utilities.LFP.defineChannelListSFC(patientIDs, cds, chanType);

%% Execute
% only collecting spike data takes ~1 min
% collecting spike and lfp data takes ~15 min
% collecting and computing with boot takes ~30min

tic   
for sess = 1:length(dirID)
    
    if sess ~= 1
        params = rmfield(params, 'valid_cell');
    end
    
    sessDir = dirID(sess, :);
    
    [lfpChans] = Utilities.LFP.defineLFPChannelsSFC(sessDir, lfpArea);
    
    % read in session IDs and assign conditions etc.
    [params] = Utilities.LFP.defineInputParamsSFC(cellArea, lfpArea, lfpChans, sessDir{1}, 'log', 'sigRamp', diskPath);
    
    params.chanType = chanType;
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
        
        % compute ppc for that session and condition
        if ~isempty(data_lfp) && ~saveDataOnly
            [ppc{sess, cond}, frq, ppc_boot{sess, cond}] = Recall.compute_ppc(data_lfp{sess, cond}, data_spike{sess, cond}, params.low_freq, params.high_freq, n_freq, params.scale, params.FsDown, params.run_boot);
        end
   
    end
    
    
    
end
toc

%% save 

if saveDataOnly
    outPath = [diskPath filesep 'Recall_task' filesep 'ppc_' params.scale filesep 'Data'];
    if ~exist(outPath, 'dir')
        mkdir(outPath)
    end
    % save data
    filename = [outPath filesep ['data_SFC_' params.cellArea 'Cell_' params.lfpArea 'LFP_' params.cells '_' cds '_' params.chanType]];
    save(filename, 'data_spike', 'data_lfp', 'params');
else
    outPath = [diskPath filesep 'Recall_task' filesep 'ppc_' params.scale ];
    if ~exist(outPath, 'dir')
        mkdir(outPath)
    end
    filename = [outPath filesep ['ppc_' params.cellArea 'Cell_' params.lfpArea 'LFP_' params.cells '_' cds  '_' params.chanType]];
    save(filename, 'ppc', 'ppc_boot', 'frq', 'params');
end
%% Go to compute_SFC_Stats or plot_SFC_main 



