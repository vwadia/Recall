% Imagined decoding
% finally getting down to this
% leggo.

% Things I will need:
% feature matrices: cells x trials per session
% labels: the ground truth image identities

% Steps I will need:
% destroy correlations between cells?
% Then decode w/ SVMDecoder (cross validated)
% or with fitcdisc (disciminant analysis)?

% Plot scatter with decoding accuracy per session (maybe across time?)
% First just across all 5s period

% Idea for testing reactivation of cells
% if ttest of baseline vs Im period is significant
    
%% Prep data from each session

% Could compute with all cells, all IT cells, all axis tuned cells

setDiskPaths

sessID = {['Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917'],...
    ['Recall_Task' filesep 'P76CS' filesep 'Recall_Session_2_20210925'],... % want the rasters so using recall session
    ['Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_3_20210927'],...
    ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_1_20220330'],...
    ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_2_20220403'],...
    ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_3_20220405'],...
    ['Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_1_20220728'],...
    ['Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_2_20220731']};



for ss = 1:length(sessID)
    
    load([diskPath filesep sessID{ss} filesep 'RecallData_NoFreeRec.mat'])
    strctCELL = struct2cell(strctCells');
    strctCELL = strctCELL';

    imRasters = RecallData.CRTimeCourse;
    
    labels = RecallData.CROrder;
    
end
