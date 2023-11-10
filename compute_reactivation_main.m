% software package versino of the ugly ass scripts 
% computeReactivation and computeReactivation _2

setDiskPaths

% grab TTLs
taskCodePath = [boxPath filesep 'RecallTaskVarun'];
addpath(taskCodePath)
setTTLCodes

% set session list
task = 'Recall_Task';
sessID = Utilities.sessionListAllTasks(task, false);
sessID{2} = ['Recall_Task' filesep 'P76CS' filesep 'Recall_Session_2_20210925']; % want the rasters so using recall session
sessID = sessID';


for ss = 1:length(sessID)
    
    sessDir = sessID{ss};
    BlineType = 1; % encoding baseline
    method = 'threshold';
    perStim = true;
    
    % input parameters
    params = Utilities.defineInputParamsReactivation(sessDir, BlineType, method, perStim);
    
    % compute reactivation
    reac_names = Utilities.computeReactivation(params);
    
    % compute enc vs im corr
    [cc_encIm, cc_scrnIm] = Utilities.computeResponseCorr(params);
    
    % compute proj vs enc/im corr
    [cc_axCorrEnc, cc_axCorrIm] = Utilities.computeAxisCorrelations(params);
end