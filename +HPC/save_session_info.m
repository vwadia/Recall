% script to save SFC indo for each session separately (both
% conditions)

% Note this should be run locally

setDiskPaths

%% Define patients, areas, and channel list 

patientIDs = {'P76CS', 'P79CS', 'P80CS'};
cellArea = 'RFFA'; lfpArea = 'RH';

condition_1 = 'Screening'; cds = ['ScreeningImagination'];
% condition_1 = 'Encoding'; cds = ['EncodingImagination'];
condition_2 = 'Imagination';

chanType = 'Cell';
% chanType = 'NoNoise';
dirID = Utilities.LFP.defineChannelListSFC(patientIDs, cds, chanType);


% save these in a loop

for sess = 1:length(dirID)
    
    sessDir = dirID(sess, :);
    
    [lfpChans] = Utilities.LFP.defineLFPChannelsSFC(sessDir, lfpArea);
    
    % read in session IDs and assign conditions etc.
    [params] = Utilities.LFP.defineInputParamsSFC(cellArea, lfpArea, lfpChans, sessDir{1}, 'log', 'sigRamp', diskPath);
    
    params.chanType = chanType; % doing this before defineInputParams doesn't work
    
    saveDir = [diskPath filesep 'Recall_Task' filesep 'SFC_Session_info'];
    if ~exist(saveDir, 'dir')
        mkdir(saveDir)
    end
    
    save([saveDir filesep 'SFC_session_' num2str(sess) '_info_HPC.mat'], 'params', 'condition_1', 'condition_2');
    
end
