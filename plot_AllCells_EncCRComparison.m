% script to concatenate all IT cells in recall (for plotting purposes)



setDiskPaths

taskCodePath = [boxPath filesep 'RecallTaskVarun'];
addpath(taskCodePath)
setTTLCodes

% set session list
task = 'Recall_Task';
sessID = Utilities.sessionListAllTasks(task, false);
sessID{2} = ['Recall_Task' filesep 'P76CS' filesep 'Recall_Session_2_20210925']; % want the rasters so using recall session

alpha = 0.05;
load([diskPath filesep task filesep ['ReactiveITCells_alpha' num2str(alpha) '_500Stim_Im_SigRamp.mat']]); 
RC = strctCells;
clearvars psths responses strctCells

% % load manually curated reactive cells
% RC_names_man = readtable([diskPath filesep task filesep 'ReactiveITCells_manual']);
% RC_nm = [table2array(RC_names_man(:, 1)), table2array(RC_names_man(:, 3))];
% RC_names = RC_nm(:, 1);
% RC_ChanNum = RC_nm(:, 2);

RC_names = cat(1, RC(:).Name);
RC_ChanNum = cat(1, RC(:).ChannelNumber);

ITOnly = true;
tic 
for ss = 1:length(sessID)
     
     
    load([sessID{ss} filesep 'RecallData_NoFreeRec.mat'])
    strctCELL = struct2cell(strctCells')';
    
    if ITOnly
        IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
        full_strctCELL = strctCELL;
        strctCELL = strctCELL(IT_Cells, :);
        strctCells = strctCells(IT_Cells);
        EncodingTimeCourse =  RecallData.EncodingTimeCourse(IT_Cells, 1:3);
        CRTimeCourse = RecallData.CRTimeCourse(IT_Cells, 1:3);
        
        % Baseline for CR
        BTimeCourse = RecallData.PreTrialBaselineTimeCourse;        
    end
    parfor cellIndex = 1:length(strctCells)
        enc_psth = EncodingTimeCourse(cellIndex, 1:3);
        CR_psth = CRTimeCourse(cellIndex, 1:3);
        EncodingOrder = RecallData.EncodingOrder;
        CROrder = RecallData.CROrder;
        offsetEnc = RecallData.offsetEnc;
        offsetTones = RecallData.offsetTones;
        stimuli = RecallData.stimuli;
        
        % subbing baseline in for FR 
%         FR_psth = BTimeCourse(cellIndex, 1:3);
%         FROrder = ones(length(RecallData.stimuli), 1);
%         offsetFR = RecallData.offsetFR;        
%         f = Utilities.Plotting.PlotPerCellAllStim_Im(enc_psth, CR_psth, EncodingOrder, CROrder, offsetEnc, offsetTones, stimuli, FR_psth, FROrder, offsetFR);
        
        f = Utilities.Plotting.PlotPerCellAllStim_Im(enc_psth, CR_psth, EncodingOrder, CROrder, offsetEnc, offsetTones, stimuli);
        
        sgt = sgtitle({[num2str(strctCells(cellIndex).Name) '\_' strctCells(cellIndex).brainArea '\_Comparison']}); % backslash allows you to print the underscore
        sgt.FontSize = 20;
        sgt.FontWeight = 'Bold';
        
        pathOut = [diskPath filesep task filesep 'EncodingandCR_comparison'];
        if ~exist(pathOut, 'dir')
            mkdir(pathOut)
        end
        
        if sum(ismember(RC_names, strctCells(cellIndex).Name)) ~= 0 && isequal(strctCells(cellIndex).ChannelNumber, RC_ChanNum(find(ismember(RC_names, strctCells(cellIndex).Name))))
            pathOut = [pathOut filesep 'ReactivatedCells'];
            if ~exist(pathOut, 'dir')
                mkdir(pathOut)
            end           
        end
        filename = [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name) '_Comparison'];
        
        print(f, filename, '-dpng','-r0');
        close all
    end
    
end
toc