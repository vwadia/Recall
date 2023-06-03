% second script to do the same for some legibility

setDiskPaths

taskCodePath = [boxPath filesep 'RecallTaskVarun'];
addpath(taskCodePath)
setTTLCodes

% set session list
task = 'Recall_Task';
sessID = Utilities.sessionListAllTasks(task, false);
sessID{2} = ['Recall_Task' filesep 'P76CS' filesep 'Recall_Session_2_20210925']; % want the rasters so using recall session
sessID = sessID';
manual_curation = false; % did I look for reactivation by eye?

if manual_curation
    % load manually curated reactive cells
    RC_names_man = readtable([diskPath filesep task filesep 'ReactiveITCells_manual.xlsx']);
    RC_nm = [table2array(RC_names_man(:, 1)), table2array(RC_names_man(:, 3))];
    
    qm_cells = [1597 2150 2927 2884 1483 2522 4732 1354 3502 5808 1569 2066 5733];
    RC_nm = RC_nm(~ismember(RC_nm(:, 1), qm_cells), :);
end


ITOnly = true;
perStim = true;
alpha = 0.05;

im_active_cells = {};
ctr = 1;
all_ctr = 1;
BlineType = 1;
perStimBase = true;
useThreshold = true;

n_stdDevs = 5;



for ss = 1:length(sessID)
    tic
    
    load([sessID{ss} filesep 'RecallData_NoFreeRec.mat'])
    
    strctCELL = struct2cell(strctCells')';
    switch BlineType
        case 1
            BTimeCourse = RecallData.EncodingTimeCourse;
            b_end = 1500;
        case 2
            BTimeCourse = RecallData.PreCRBaselineTimeCourse;
            b_end = 5000;
        case 3
            BTimeCourse = RecallData.PreTrialBaselineTimeCourse;
            b_end = 5000;
    end
    
    CRTimeCourse = RecallData.CRTimeCourse;
    
    if ITOnly
        IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
        full_strctCELL = strctCELL;
        strctCELL = strctCELL(IT_Cells, :);
        strctCells = strctCells(IT_Cells);
        
        BTimeCourse = BTimeCourse(IT_Cells, :);
        CRTimeCourse = CRTimeCourse(IT_Cells, :);
    end
    
    sigcells = false(1, length(strctCells));
    
    % find max group and then check whether it exceeds baseline
    for cellIndex = l(strctCells)
        
        % grab baseline
        if ~perStimBase || BlineType ~= 1
            if BlineType == 1
                b_psth = BTimeCourse{cellIndex, 2}(:, 1:b_end);                
            else
                b_psth = BTimeCourse{cellIndex, 2}(:, 2000:b_end);
            end
            b_timecourse = mean(b_psth, 1);
            h_threshold = mean(b_timecourse) + (n_stdDevs*std(b_timecourse));
            l_threshold = mean(b_timecourse) + (3*std(b_timecourse));
        end
        CROrder = RecallData.CROrder;
        max_cr = nan(1, length(unique(CROrder)));
        p = nan(1, length(unique(CROrder')));
        for stim = unique(CROrder')
            if perStimBase && BlineType == 1
                b_psth = BTimeCourse{cellIndex, 2}(CROrder == stim, 1:b_end);
                b_timecourse = mean(b_psth, 1);
                h_threshold = mean(b_timecourse) + (n_stdDevs*std(b_timecourse));
                l_threshold = mean(b_timecourse) + (3*std(b_timecourse));
                
                b_trials = mean(BTimeCourse{cellIndex, 1}(CROrder == stim, 1:b_end), 2);
                
            end
            cr_psth = CRTimeCourse{cellIndex, 2}(CROrder == stim, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000);
            cr_timecourse = mean(cr_psth, 1);
            max_cr(stim) = max(cr_timecourse);
            
            cr_trials = mean(CRTimeCourse{cellIndex, 1}...
                (CROrder == stim, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000), 2);
            
            if exist('b_trials', 'var')
                [p(stim), ~] = ranksum(b_trials, cr_trials, 'tail', 'left');
            end
         
        end
        if useThreshold
            [max_val, max_stim] = max(max_cr);
            
            if max_val > h_threshold %|| (sum(max_cr > l_threshold) == length(unique(CROrder)))
                %         if (sum(max_cr > l_threshold) == length(unique(CROrder)))
                sigcells(cellIndex) = true;
                im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
            end
        else
            if sum(p < alpha/length(unique(CROrder'))) > 0
                sigcells(cellIndex) = true;
                im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
            end
        end
    end
        
           
    fprintf('Finished for session %d', ss)
    fprintf('\n')
    toc
       
        
end
 
    
%% compare to axis tuning

load([diskPath filesep task filesep 'AllITCells_500Stim_Im_SigRamp.mat'])

if ~manual_curation
    reac_names = [cell2mat(im_active_cells(:, 1)) cell2mat(im_active_cells(:, 2))];
else
    reac_names = RC_nm;
end
sigramp_names = [cat(1, strctCells(:).Name), cat(1, strctCells(:).ChannelNumber)];

% which cells did both
ovrlap = ismember(sigramp_names(:, 1), reac_names(:, 1));
% make sure channel numbers are the same too
for i = 1:length(ovrlap)
    
    if ovrlap(i) == 1
        rel_idx = find(reac_names(:, 1) == sigramp_names(i, 1));
        
        if ~isequal(sigramp_names(i, 2), reac_names(rel_idx, 2))
            ovrlap(i) = false;
        end
    end
    
end

strctCells = strctCells(ovrlap);
responses = responses(ovrlap, :);
psths = psths(ovrlap, :);

% save([diskPath filesep task filesep ['SigRampCellsthatReactivate_alpha' num2str(alpha) '.mat']], 'ovrlap')
% save([diskPath filesep task filesep ['ReactiveITCells_alpha ' num2str(alpha) '_500Stim_Im_SigRamp.mat']], 'strctCells', 'responses', 'psths');


%%

% which cells did both
check_ovrlap = ismember(reac_names(:, 1), RC_nm(:, 1));
% make sure channel numbers are the same too
for i = 1:length(check_ovrlap)
    
    if check_ovrlap(i) == 1
        rel_idx = find(reac_names(i, 1) == RC_nm(:, 1));
        
        if ~isequal(RC_nm(rel_idx, 2), reac_names(i, 2))
            check_ovrlap(i) = false;
        end
    end
    
end

chosenCells = im_active_cells(check_ovrlap, :);
nonchosenCells = im_active_cells(~check_ovrlap, :);
% responses = responses(ovrlap, :);
% psths = psths(ovrlap, :);



