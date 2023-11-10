% script to compute reactivation metric

% best way to do this?
% loop through all sessions
% per session
% per cell
%    Baseline options
%       Encoding baseline (DONE)
%       Time before Imagination period (Need to remake psths to accommodate this)
%       Period between trials (Need to remake psths to accommodate this)
% plot rasters to include that baseline

% select cells that activated in Im
% choose overlap of cells that had vis resp and reactivated in Im

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
    RC_names_man = readtable([diskPath filesep task filesep 'ReactiveITCells_manual']);
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


useCombo = false;
useAnova = true;
useThreshold = false;
n_stdDevs = 5;
full_strctCELL = {};

includeScreening  = false;


for ss = 1:length(sessID)
    tic
    %     if ss ~= 2
    load([sessID{ss} filesep 'RecallData_NoFreeRec.mat'])
    if includeScreening
        load([sessID{ss} filesep 'PsthandResponses.mat']);
    end
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
    EncTimeCourse = RecallData.EncodingTimeCourse;
    if ITOnly
        IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
        strctCELL = strctCELL(IT_Cells, :);
        strctCells = strctCells(IT_Cells);
        
        BTimeCourse = BTimeCourse(IT_Cells, :);
        CRTimeCourse = CRTimeCourse(IT_Cells, :);
        
        if exist('responses', 'var')
            responses = responses(IT_Cells, :);
            psths = psths(IT_Cells, :);
            labels_train = order;
        end
        
        EncTimeCourse = RecallData.EncodingTimeCourse(IT_Cells, :);
    end
    
    
    if exist('labels_train', 'var')
        % covert train and test labels to match
        train_files = Utilities.readInFiles([diskPath filesep sessID{ss} filesep 'stimuliUsed'], 'tif');
        train_files = struct2cell(train_files)';
        test_files = Utilities.readInFiles([diskPath filesep sessID{ss} filesep 'stimuliUsedRecall'], 'tif');
        test_files = struct2cell(test_files)';
        train_files = cellfun(@(x) x(1:end-4), train_files(:, 1), 'UniformOutput', false);
        test_files = cellfun(@(x) x(1:end-4), test_files(:, 1), 'UniformOutput', false);
        
        matched = find(ismember(train_files, test_files) == 1);
        
    end
    
    full_strctCELL = cat(1, full_strctCELL, strctCELL);
    
    % concatenate anova and others
    % does best 37/96, 0.39 and good corr (with bline 3)
    if useCombo
        
        num_cells(ss) = length(strctCells);
        sigcells = false(1, length(strctCells));
        for cellIndex = l(strctCells)
            cr_psth = CRTimeCourse{cellIndex, 1}(:, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000);
            cr_trials = mean(cr_psth, 2);
            
            group = RecallData.CROrder;
            p = anova1(cr_trials, group, 'off');
            
            if p < alpha
                sigcells(cellIndex) = true;
                im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
            end
            
            % for correlation computation
            if ss == 1
                idx = cellIndex;
            else
                idx = sum(num_cells(1:ss-1)) + cellIndex;
            end
            
            % grab the encoding responses - like in plotPerCellRasterandbarPerStim
            if exist('responses', 'var') && ~isnan(responses{cellIndex, 2})
                start_e = RecallData.offsetEnc(1) + floor(responses{cellIndex, 2});
                till_e = start_e + 500;%RecallData.stimDur*1e3;
                
                start_s = floor(responses{cellIndex, 2}) + 170;
                till_s = start_s + 267;
                
            else
                start_e = RecallData.offsetEnc(1);
                till_e = start_e + 500;% RecallData.stimDur*1e3;
                
                start_s = 170;
                till_s = start_s + 267;
                
            end
            start_cr = RecallData.offsetTones(1);
            till_cr = RecallData.offsetTones(1)+5000;
            for stim = unique(RecallData.EncodingOrder')
                full_EncResponses{idx, 1}(stim, 1) = mean(mean(EncTimeCourse{cellIndex, 1}(find(RecallData.EncodingOrder == stim), start_e:till_e)))*1e3;
                full_EncResponses{idx, 2} = strctCells(cellIndex).Name;
            end
            if includeScreening
                ctr_stim = 1;
                for stim = unique(matched')
                    full_ScrnResponses{idx, 1}(ctr_stim, 1) = mean(mean(psths{cellIndex, 1}(find(order == stim), start_s:till_s)))*1e3;
                    ctr_stim = ctr_stim+1;
                end
            end
            for stim = unique(RecallData.CROrder')
                full_CRResponses{idx, 1}(stim, 1) = mean(mean(CRTimeCourse{cellIndex, 1}(find(RecallData.CROrder == stim), start_cr:till_cr)))*1e3;
            end
            
        end
        
        
        % take other neurons
        strctCells = strctCells(~sigcells);
        strctCELL = strctCELL(~sigcells, :);
        BTimeCourse = BTimeCourse(~sigcells, :);
        CRTimeCourse = CRTimeCourse(~sigcells, :);
        
        for cellIndex = l(strctCells)
            if useThreshold
                threshCount = 0;
                
                b_psth = BTimeCourse{cellIndex, 2}(:, 1:b_end);
                b_timecourse = mean(b_psth, 1);
                threshold = mean(b_timecourse) + (n_stdDevs*std(b_timecourse));
                if perStim
                    CROrder = RecallData.CROrder;
                    p = nan(1, length(unique(CROrder)));
                    max_cr = nan(1, length(unique(CROrder)));
                    for stim = unique(CROrder')
                        cr_psth = CRTimeCourse{cellIndex, 2}(CROrder == stim, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000);
                        cr_timecourse = mean(cr_psth, 1);
                        max_cr(stim) = max(cr_timecourse);
                        if max(cr_timecourse) > threshold
                            %                             threshCount = threshCount + 1;
                            
                            im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
                            break;
                            
                        end
                    end
                    %                     if threshCount >= length(unique(CROrder'))/2 || max(max_cr) >  mean(b_timecourse) + (5*std(b_timecourse))
                    %                         im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
                    %                     end
                else
                    cr_psth = CRTimeCourse{cellIndex, 2}(:, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000);
                    cr_timecourse = mean(cr_psth, 1);
                    if max(cr_timecourse) > threshold
                        im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
                    end
                end
            else
                % get cell baseline during encoding
                b_psth = BTimeCourse{cellIndex, 1}(:, 1:b_end);
                
                b_trials = mean(b_psth, 2);
                
                % grab the imagination rasters *per stim*
                if perStim
                    CROrder = RecallData.CROrder;
                    p = nan(1, length(unique(CROrder)));
                    for stim = unique(CROrder')
                        cr_psth = CRTimeCourse{cellIndex, 1}(CROrder == stim, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000);
                        cr_trials = mean(cr_psth, 2);
                        [p(stim), ~] = ranksum(b_trials, cr_trials, 'tail', 'left'); % checking that cr_trials is bigger than b_trials not just different
                    end
                    if sum(p < alpha/length(p)) ~= 0 % bonferroni correction
                        %                 if sum(p < alpha) ~= 0
                        im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
                    end
                    
                end
            end
            
            
            %             full_CRTimeCourse = cat(1, full_CRTimeCourse, CRTimeCourse(cellIndex, :));
            %             full_EncTimeCourse = cat(1, full_EncTimeCourse, RecallData.EncodingTimeCourse(cellIndex, :)); % grab from reacllData so you can use it with diff baselines
        end
        
        
    else
        
        % for each cell, grab baseline across all stim
        % compare to activity in Im across all stim
        % Could also compare to activity in Im per stim (keep cells that have
        % sig act for at  least 1 stim)
        num_cells(ss) = length(strctCells);
        for cellIndex = l(strctCells)
            
            
            if ~useAnova
                
                if useThreshold
                    threshCount = 0;
                    
                    b_psth = BTimeCourse{cellIndex, 2}(:, 2000:b_end);
                    b_timecourse = mean(b_psth, 1);
                    threshold = mean(b_timecourse) + (n_stdDevs*std(b_timecourse));
                    if perStim
                        CROrder = RecallData.CROrder;
                        p = nan(1, length(unique(CROrder)));
                        max_cr = nan(1, length(unique(CROrder)));
                        for stim = unique(CROrder')
                            cr_psth = CRTimeCourse{cellIndex, 2}(CROrder == stim, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000);
                            cr_timecourse = mean(cr_psth, 1);
                            max_cr(stim) = max(cr_timecourse);
                            if max(cr_timecourse) > threshold
                                threshCount = threshCount + 1;
                                %                             if threshCount >= length(unique(CROrder'))/2
                                %                                 im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
                                %                                 break;
                                %                             end
                            end
                        end
                        if threshCount >= length(unique(CROrder'))/2 || max(max_cr) >  mean(b_timecourse) + (5*std(b_timecourse))
                            im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
                        end
                    else
                        cr_psth = CRTimeCourse{cellIndex, 2}(:, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000);
                        cr_timecourse = mean(cr_psth, 1);
                        if max(cr_timecourse) > threshold
                            im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
                        end
                    end
                else
                    % get cell baseline during encoding
                    b_psth = BTimeCourse{cellIndex, 1}(:, 1:b_end);
                    
                    b_trials = mean(b_psth, 2);
                    
                    % grab the imagination rasters *per stim*
                    if perStim
                        CROrder = RecallData.CROrder;
                        p = nan(1, length(unique(CROrder)));
                        for stim = unique(CROrder')
                            cr_psth = CRTimeCourse{cellIndex, 1}(CROrder == stim, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000);
                            cr_trials = mean(cr_psth, 2);
                            [p(stim), ~] = ranksum(b_trials, cr_trials, 'tail', 'left'); % checking that cr_trials is bigger than b_trials not just different
                        end
                        
                        all_pvals{all_ctr, :} = p;
                        all_ctr = all_ctr+1;
                        
                        if sum(p < alpha/length(p)) ~= 0 % bonferroni correction
                            %                 if sum(p < alpha) ~= 0
                            im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
                            pvals{ctr, :} = p;
                            ctr = ctr + 1;
                        end
                        
                    else
                        cr_psth = CRTimeCourse{cellIndex, 1}(:, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000);
                        cr_trials = mean(cr_psth, 2);
                        
                        [p, h] = ranksum(b_trials, cr_trials, 'tail', 'left');
                        %         [h, p] = ttest(b_trials, cr_trials);
                        
                        all_pvals(all_ctr, 1) = p;
                        all_ctr = all_ctr+1;
                        if p < alpha
                            im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
                            pvals(ctr, 1) = p;
                            ctr = ctr + 1;
                        end
                    end
                end
            else
                
                cr_psth = CRTimeCourse{cellIndex, 1}(:, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000);
                cr_trials = mean(cr_psth, 2);
                
                group = RecallData.CROrder;
                p = anova1(cr_trials, group, 'off');
                
                all_pvals(all_ctr, 1) = p;
                all_ctr = all_ctr+1;
                if p < alpha
                    im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
                    pvals(ctr, 1) = p;
                    ctr = ctr + 1;
                end
                
                
                % for correlation computation
                if ss == 1
                    idx = cellIndex;
                else
                    idx = sum(num_cells(1:ss-1)) + cellIndex;
                end
                % grab the encoding responses - like in plotPerCellRasterandbarPerStim
                if exist('responses', 'var') && ~isnan(responses{cellIndex, 2})
                    start_e = RecallData.offsetEnc(1) + floor(responses{cellIndex, 2});
                    till_e = start_e + 500;%RecallData.stimDur*1e3;
                    
                    start_s = floor(responses{cellIndex, 2}) + 170;
                    till_s = start_s + 267;
                    
                else
                    start_e = RecallData.offsetEnc(1);
                    till_e = start_e + 500;% RecallData.stimDur*1e3;
                    
                    start_s = 170;
                    till_s = start_s + 267;
                    
                end
                start_cr = RecallData.offsetTones(1);
                till_cr = RecallData.offsetTones(1)+5000;
                for stim = unique(RecallData.EncodingOrder')
                    full_EncResponses{idx, 1}(stim, 1) = mean(mean(EncTimeCourse{cellIndex, 1}(find(RecallData.EncodingOrder == stim), start_e:till_e)))*1e3;
                    full_EncResponses{idx, 2} = strctCells(cellIndex).Name;
                end
                if includeScreening
                    ctr_stim = 1;
                    for stim = unique(matched')
                        full_ScrnResponses{idx, 1}(ctr_stim, 1) = mean(mean(psths{cellIndex, 1}(find(order == stim), start_s:till_s)))*1e3;
                        ctr_stim = ctr_stim+1;
                    end
                end
                for stim = unique(RecallData.CROrder')
                    full_CRResponses{idx, 1}(stim, 1) = mean(mean(CRTimeCourse{cellIndex, 1}(find(RecallData.CROrder == stim), start_cr:till_cr)))*1e3;
                end
            end
        end
        
    end
    
    fprintf('Finished for session %d', ss)
    fprintf('\n')
    toc
    %     end
    
end
full_strctCells = cell2struct(full_strctCELL, fieldnames(strctCells), 2);

% if iscell(pvals)
%     min_p = cell2mat(cellfun(@(x) min(x), pvals(:, 1), 'UniformOutput', false));
% end
% tt1 = cellfun(@(x) mean(x)*1e-3, times(:, 1), 'UniformOutput', false);
% tt2 = cellfun(@(x) mean(x)*1e-3, times(:, 2), 'UniformOutput', false);

%% load in axis tuned neurons

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

% strctCells = strctCells(ovrlap);
% responses = responses(ovrlap, :);
% psths = psths(ovrlap, :);
%
% if useCombo save([diskPath filesep task filesep ['SigRampCellsthatReactivate_alpha' num2str(alpha) '_Combo.mat']], 'ovrlap');
%     save([diskPath filesep task filesep ['ReactiveITCells_alpha' num2str(alpha) '_500Stim_Im_SigRamp_Combo.mat']], 'strctCells', 'responses', 'psths');
%
%
% else
%     save([diskPath filesep task filesep ['SigRampCellsthatReactivate_alpha' num2str(alpha) '.mat']], 'ovrlap');
%     save([diskPath filesep task filesep ['ReactiveITCells_alpha' num2str(alpha) '_500Stim_Im_SigRamp.mat']], 'strctCells', 'responses', 'psths');
%
% end

%% pulling out various combinations of cells
full_names = [cell2mat(full_strctCELL(:, 1)) cell2mat(full_strctCELL(:, 2))];

% axis tuned neurons
load([diskPath filesep task filesep 'AllITCells_500Stim_Im_SigRamp.mat'])
sigramp_names = [cat(1, strctCells(:).Name), cat(1, strctCells(:).ChannelNumber)];

% which cells did both
ovrlap_ax = ismember(full_names(:, 1), sigramp_names(:, 1));
% make sure channel numbers are the same too
for i = 1:length(ovrlap_ax)
    
    if ovrlap_ax(i) == 1
        rel_idx = find(sigramp_names(:, 1) == full_names(i, 1));
        
        if ~isequal(full_names(i, 2), sigramp_names(rel_idx, 2))
            ovrlap_ax(i) = false;
        end
    end
    
end



% responsive neurons
load([diskPath filesep task filesep 'AllRespITCells_500Stim_Im.mat'])
R_strctCELL = struct2cell(strctCells')';
resp_names = [cell2mat(R_strctCELL(:, 1)) cell2mat(R_strctCELL(:, 2))];

% which cells did both
ovrlap_resp = ismember(full_names(:, 1), resp_names(:, 1));
% make sure channel numbers are the same too
for i = 1:length(ovrlap_resp)
    
    if ovrlap_resp(i) == 1
        rel_idx = find(resp_names(:, 1) == full_names(i, 1));
        
        if ~isequal(full_names(i, 2), resp_names(rel_idx, 2))
            ovrlap_resp(i) = false;
        end
    end
    
end

% reactive neurons
if ~manual_curation
    reac_names = [cell2mat(im_active_cells(:, 1)) cell2mat(im_active_cells(:, 2))];
else
    reac_names = RC_nm;
end
% which cells did both
ovrlap = ismember(full_names(:, 1), reac_names(:, 1));
% make sure channel numbers are the same too
for i = 1:length(ovrlap)
    
    if ovrlap(i) == 1
        rel_idx = find(reac_names(:, 1) == full_names(i, 1));
        
        if ~isequal(full_names(i, 2), reac_names(rel_idx, 2))
            ovrlap(i) = false;
        end
    end
    
end

if includeScreening
    reacData = table(full_names(:, 1), full_names(:, 2), cc_all', cc_all_s', ovrlap, ovrlap_resp, ovrlap_ax);
    reacData.Properties.VariableNames = {'Name', 'Channel', 'encCorr','scrnCorr','reactivated', 'responsive', 'axisTuned'};
    % save([diskPath filesep task filesep 'ReactivationCorrelationData_PerStimBaseThresh_wScrn' num2str(0.05) '.mat'], 'reacData');
else
    reacData = table(full_names(:, 1), full_names(:, 2), cc_all', ovrlap, ovrlap_resp, ovrlap_ax);
    reacData.Properties.VariableNames = {'Name', 'Channel', 'encCorr','reactivated', 'responsive', 'axisTuned'};
    % save([diskPath filesep task filesep 'ReactivationCorrelationData_PerStimBaseThresh' num2str(0.05) '.mat'], 'reacData');
    
end
reacData = sortrows(reacData, 'encCorr', 'descend');
%%
% load([diskPath filesep task filesep 'AllITCells_500Stim_Im_SigRamp.mat'])
% sigramp_names = [cat(1, strctCells(:).Name), cat(1, strctCells(:).ChannelNumber)];
%
% % check for names that overlap
% ovrlap = ismember(sigramp_names(:, 1), RC_nm(:, 1));
%
% % make sure channel numbers are the same too
% for i = 1:length(ovrlap)
%
%     if ovrlap(i) == 1
%         rel_idx = find(RC_nm(:, 1) == sigramp_names(i, 1));
%
%         if ~isequal(sigramp_names(i, 2), RC_nm(rel_idx, 2))
%             ovrlap(i) = false;
%         end
%     end
%
% end






