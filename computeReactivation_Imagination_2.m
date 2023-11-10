% second script to do the same for some legibility
% use this one for my agreed upon imagination metric
% vwadia June2023
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
incAnova = true;

n_stdDevs = 5;

includeScreening  = true;

% for correlation computation
full_CRTimeCourse = {};
full_EncTimeCourse = {};
full_strctCELL = {};
full_EncResponses = {};
full_CRResponses = {};
full_strctCells = struct;

if includeScreening
    full_ScrnResponses = {};
    full_ITResp = {};
end


for ss = 1:length(sessID) 
    tic

    if includeScreening && ss == 2 % this session didn't have screening
        continue
    end 
    load([sessID{ss} filesep 'RecallData_NoFreeRec.mat'])
    if includeScreening 
        load([sessID{ss} filesep 'PsthandResponses.mat']);        
        load([sessID{ss} filesep 'ITResponses.mat']) % strctResp
        full_ITResp = cat(1, full_ITResp, struct2cell(strctResp')');
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
%         full_strctCELL = strctCELL;
        strctCELL = strctCELL(IT_Cells, :);
        strctCells = strctCells(IT_Cells);
        
        BTimeCourse = BTimeCourse(IT_Cells, :);
        CRTimeCourse = CRTimeCourse(IT_Cells, :);
        
        EncTimeCourse = RecallData.EncodingTimeCourse(IT_Cells, :);
        
        if exist('responses', 'var') && ss ~= 2
            responses = responses(IT_Cells, :);
            psths = psths(IT_Cells, :);
            labels_train = order;
        end
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
    
    sigcells = false(1, length(strctCells));
    num_cells(ss) = length(strctCells);

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
            l_threshold = mean(b_timecourse) + (2*std(b_timecourse));
        end
        
        CROrder = RecallData.CROrder;
        EncOrder = RecallData.EncodingOrder;
        max_cr = nan(1, length(unique(CROrder)));
        p = nan(1, length(unique(CROrder')));
        if exist('b_psth', 'var')
            b_trials = mean(b_psth, 2);
        end
        
        full_b_timecourse = mean(BTimeCourse{cellIndex, 2}(:, 1:b_end), 1);
        full_h_threshold = mean(full_b_timecourse) + n_stdDevs*std(full_b_timecourse);
        full_l_threshold = mean(full_b_timecourse) + 2*std(full_b_timecourse);
        full_b_trials = mean(BTimeCourse{cellIndex, 1}(:, 1:b_end), 2);
        
        if incAnova
            cr_psth = CRTimeCourse{cellIndex, 1}(:, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000);
            cr_trials = mean(cr_psth, 2);
            
            group = RecallData.CROrder;
            p_anov = anova1(cr_trials, group, 'off');
 
        end
        
        % compute values 
        for stim = unique(CROrder')
            
            e_psth = EncTimeCourse{cellIndex, 2}(EncOrder == stim, RecallData.offsetEnc(1):end);
            e_timecourse = mean(e_psth, 1);
            e_trials = mean(EncTimeCourse{cellIndex, 1}(EncOrder == stim, RecallData.offsetEnc(1):end), 2);
            
            if perStimBase && BlineType == 1
                b_psth = BTimeCourse{cellIndex, 2}(CROrder == stim, 1:b_end);
                b_timecourse = mean(b_psth, 1);
                h_threshold(stim) = mean(b_timecourse) + (n_stdDevs*std(b_timecourse));
                l_threshold(stim) = mean(b_timecourse) + (2*std(b_timecourse));
                b_mean(stim) = mean(b_timecourse);
                b_trials = mean(BTimeCourse{cellIndex, 1}(CROrder == stim, 1:b_end), 2);
           
            end
            cr_psth = CRTimeCourse{cellIndex, 2}(CROrder == stim, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000);
            cr_timecourse = mean(cr_psth, 1);
            max_cr(stim) = max(cr_timecourse);
            mean_cr(stim) = mean(cr_timecourse);
            
            cr_trials = mean(CRTimeCourse{cellIndex, 1}...
                (CROrder == stim, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000), 2);
            
            if exist('b_trials', 'var')
                [p_2(stim), ~] = ranksum(b_trials, cr_trials, 'tail', 'left');
                [p(stim), ~] = ranksum(full_b_trials, cr_trials, 'tail', 'left');
            end
         
        end
        

        
        % check for reac
        if useThreshold
            [max_val, max_stim] = max(max_cr);
            
%             if sum(max_cr > mean(full_b_timecourse)) > 0
%            if sum(max_cr > full_h_threshold) > 0  % 81/96 
%            if max_val > full_h_threshold % way too permissive
%             if max_val > h_threshold(max_stim) % 47/96 but terrible ax corr 0.21 same as all cells
%             if sum(mean_cr > mean(b_timecourse)) > 0 % 82/96 could work as a pre-screen
%             if sum(mean_cr > mean(full_b_timecourse)) > 0 % 83/96 could work as a pre-screen                
%             if (sum(mean_cr > l_threshold) > 0) || max_val > h_threshold(max_stim) % pretty permissive 54/96, but ax corr 0.24 and good ecorr
%             if (sum(mean_cr > full_l_threshold) > 0) || max_val > h_threshold(max_stim) % 47/96 not great ax corr 0.21 (same as all) but good ecorr
%             if (sum(max_cr > h_threshold) > 0)  % way too permissive
%             if (sum(mean_cr > l_threshold) > 0) % 25/96, decent ecorr and great ax corr 0.3
%             if (sum(mean_cr > full_l_threshold) > 1) ||  max_val > h_threshold(max_stim) % 46/96, misses key cells and awful ax corr 0.20
%             if (sum(mean_cr > full_l_threshold) > 0) ||  max_val > h_threshold(max_stim) % 65/96, ok ax corr 0.25
%             if (sum(mean_cr > full_l_threshold) > 0) % 19/96, decent ecorr, good ax corr 0.27
%             if (sum(mean_cr > full_l_threshold) > 0) || p_anov < alpha % 39/96, decent ecorr, great ax corr 0.3 - good option but(non-reac cells have 0.16 ax corr - sig)
%             if (sum(mean_cr > full_l_threshold) > 0) && p_anov < alpha % 6/96 LOL but amaze ax corr 0.56
            if (sum(mean_cr > l_threshold) > 0) || p_anov < alpha % 44/96, decent ecorr, great ax corr 0.33 - good option (non-reac cells have 0.12 ax corr - sig)
%             if (sum(mean_cr > l_threshold) > 0) && p_anov < alpha % 7/96 and awful ax corr 0.19
%             if sum(mean_cr > mean(b_timecourse)) > 0 && p_anov < alpha % 24/96 and great ax corr 0.38 - good option but(non-reac cells have 0.16 ax corr - sig)
%             if sum(mean_cr > mean(full_b_timecourse)) > 0 && p_anov < alpha % 24/96 (actually the same as above)and great ax corr 0.38
%             if (sum(mean_cr > l_threshold) > 0) || sum(p < alpha/length(unique(CROrder'))) > 0 % 32/96 ok ecorr but good ax corr 0.31 - good option but(non-reac cells have 0.16 ax corr - sig)
%             if (sum(max_cr > h_threshold) > 0) || sum(p < alpha/length(unique(CROrder'))) > 0 % 74/96 
%             if max_val > h_threshold && p(max_stim) < alpha/length(unique(CROrder')) % 
%             if max_val > h_threshold(max_stim) && mean_cr(max_stim)  > l_threshold(max_stim) % 9/96
%             if max_val > full_h_threshold && mean_cr(max_stim) > full_l_threshold % 15/96 decent ax corr 0.24

                sigcells(cellIndex) = true;
                im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
%                 full_CRTimeCourse = cat(1, full_CRTimeCourse, CRTimeCourse(cellIndex, :));
%                 full_EncTimeCourse = cat(1, full_EncTimeCourse, RecallData.EncodingTimeCourse(cellIndex, :)); % grab from reacllData so you can use it with diff baselines
            end
        else
%             if sum(p < alpha/length(unique(CROrder'))) > 0 || sum(p_2 < alpha/length(unique(CROrder'))) > 0 % 32/92 low ax corr 0.23
%             if sum(p_2 < alpha/length(unique(CROrder'))) > 0 || p_anov < alpha % 45/96 ok ax corr 0.26
%             if sum(p < alpha/length(unique(CROrder'))) > 0 || p_anov < alpha % 39/96 great ax corr 0.3
%            if sum(p < alpha/length(unique(CROrder'))) > 0 % 19/92 ok ax corr 0.26
%             if sum(p < alpha/length(unique(CROrder'))) > 0 && p_anov < alpha % 6/96 amaze ax corr 0.56
            if sum(p_2 < alpha/length(unique(CROrder'))) > 0 && p_anov < alpha % 9/96 ok amaze ax corr 0.45

                sigcells(cellIndex) = true;
                im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
%                 full_CRTimeCourse = cat(1, full_CRTimeCourse, CRTimeCourse(cellIndex, :));
%                 full_EncTimeCourse = cat(1, full_EncTimeCourse, RecallData.EncodingTimeCourse(cellIndex, :)); % grab from reacllData so you can use it with diff baselines
            end
        end
        
        % for correlation computation
        if ss == 1
            idx = cellIndex;
        else
            idx = sum(num_cells(1:ss-1)) + cellIndex;
        end
      
       % grab the encoding responses - like in plotPerCellRasterandbarPerStim
       if exist('responses', 'var') && ~isnan(responses{cellIndex, 2}) && ss ~=2
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
       if includeScreening && ss ~= 2
           ctr_stim = 1;
           for stim = unique(matched')
               full_ScrnResponses{idx, 1}(ctr_stim, 1) = mean(mean(psths{cellIndex, 1}(find(order == stim), start_s:till_s)))*1e3;
               ctr_stim = ctr_stim+1;
           end
       end
       for stim = unique(RecallData.CROrder')
           full_CRResponses{idx, 1}(stim, 1) = mean(mean(CRTimeCourse{cellIndex, 1}(find(RecallData.CROrder == stim), start_cr:till_cr)))*1e3;
       end
         %}
    end
        
    full_strctCELL = cat(1, full_strctCELL, strctCELL);
    fprintf('Finished for session %d', ss)
    fprintf('\n')
    toc
    
end
full_strctCells = cell2struct(full_strctCELL, fieldnames(strctCells), 2);


fprintf('BlineType: %d \n', BlineType)
fprintf('perStimBase: %d \n', perStimBase)
fprintf('useThreshold: %d \n', useThreshold)
fprintf('Anova Included?: %d \n', incAnova)
fprintf('Number of Reac Cells: %d \n', length(im_active_cells))


%% compute significance of correlation - per cell
n_reps = 1000;
tic
for cellIndex = 1:length(full_strctCells)
    
    enc_s = full_EncResponses{cellIndex, 1};
    CR_s = full_CRResponses{cellIndex, 1};
    
    
    % compute correlation
    cc_boot = zeros(n_reps, 1);
    for itr = 1:n_reps
        cc_boot(itr, 1) = corr(Utilities.Shuffle(enc_s), CR_s);
%         cc_boot(itr, 1) = corr(Utilities.Shuffle(enc_s), CR_s, 'Type', 'Spearman');
    end
    
    cc_all(cellIndex) = corr(enc_s, CR_s);
%     cc_all(cellIndex) = corr(enc_s, CR_s, 'Type', 'Spearman');
    p_cc(cellIndex) = sum(cc_boot > cc_all(cellIndex))/length(cc_boot);
    
    if includeScreening
        scrn_s = full_ScrnResponses{cellIndex, 1};
        cc_boot_s = zeros(n_reps, 1);
        for itr = 1:n_reps
%             cc_boot_s(itr, 1) = corr(Utilities.Shuffle(scrn_s), CR_s, 'Type', 'Spearman');
            cc_boot_s(itr, 1) = corr(Utilities.Shuffle(scrn_s), CR_s);
        end
        cc_all_s(cellIndex) = corr(scrn_s, CR_s);
%         cc_all_s(cellIndex) = corr(scrn_s, CR_s, 'Type', 'Spearman');
        p_cc_s(cellIndex) = sum(cc_boot_s > cc_all_s(cellIndex))/length(cc_boot_s);
    end
    
    
end

fn = fieldnames(strctCells);
RC = cell2struct(im_active_cells, fn(1:6), 2);
toc
%% compare reactivation to axis tuning

% load([diskPath filesep task filesep 'AllITCells_500Stim_Im_SigRamp.mat'])
% 
% if ~manual_curation
%     reac_names = [cell2mat(im_active_cells(:, 1)) cell2mat(im_active_cells(:, 2))];
% else
%     reac_names = RC_nm;
% end
% sigramp_names = [cat(1, strctCells(:).Name), cat(1, strctCells(:).ChannelNumber)];
% 
% % which cells did both
% ovrlap = ismember(sigramp_names(:, 1), reac_names(:, 1));
% % make sure channel numbers are the same too
% for i = 1:length(ovrlap)
%     
%     if ovrlap(i) == 1
%         rel_idx = find(reac_names(:, 1) == sigramp_names(i, 1));
%         
%         if ~isequal(sigramp_names(i, 2), reac_names(rel_idx, 2))
%             ovrlap(i) = false;
%         end
%     end
%     
% end
% 
% keyboard


% strctCells = strctCells(ovrlap);
% responses = responses(ovrlap, :);
% psths = psths(ovrlap, :);

% save([diskPath filesep task filesep ['SigRampCellsthatReactivate_alpha' num2str(alpha) '.mat']], 'ovrlap')
% save([diskPath filesep task filesep ['ReactiveITCells_alpha ' num2str(alpha) '_500Stim_Im_SigRamp.mat']], 'strctCells', 'responses', 'psths');

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

fname = ['ReactivationCorrelationData' num2str(alpha)];
if includeScreening
    reacData = table(full_names(:, 1), full_names(:, 2), full_strctCELL(:, 4), cc_all', cc_all_s', ovrlap, ovrlap_resp, ovrlap_ax);
    reacData.Properties.VariableNames = {'Name', 'Channel', 'Area', 'encCorr','scrnCorr','reactivated', 'responsive', 'axisTuned'};
    fname = [fname '_wScrn'];
else
    reacData = table(full_names(:, 1), full_names(:, 2), full_strctCELL(:, 4), cc_all', ovrlap, ovrlap_resp, ovrlap_ax);
    reacData.Properties.VariableNames = {'Name', 'Channel', 'Area','encCorr','reactivated', 'responsive', 'axisTuned'};   
end
if ITOnly
    fname = [fname '_ITOnly'];
end
% reacData = sortrows(reacData, 'encCorr', 'descend');

% save([diskPath filesep task filesep fname '.mat'], 'reacData');

%%
% load([diskPath filesep task filesep 'ReactivationCorrelationData_PerStimBaseThresh' num2str(0.05) '.mat'])
% ovrlap = reacData.reactivated;

%% finding distributions for table

% % regs = unique(reacData.Area);
% regs = unique(im_active_cells(:, 4));
% 
% 
% for i = 1:numel(regs)
%     ar_name = regs{i};
%     ar_num = sum(ismember(im_active_cells(:, 4), regs{i}));
%     fprintf([ar_name ' ' num2str(ar_num) '\n'])
% 
% end



%% Corr plots

f = figure; clf;
hold on
IT_EC = reacData.encCorr(reacData.reactivated & reacData.axisTuned);
IT_EC = reacData.encCorr(reacData.reactivated & reacData.axisTuned);
histogram(IT_EC, 8, 'FaceColor', [0.4940 0.1840 0.5560]); stdErr = std(IT_EC)/sqrt(length(IT_EC));
% histogram(IT_EC, 'FaceColor', [0.4940 0.1840 0.5560]); stdErr = std(IT_EC)/sqrt(length(IT_EC));

str = {'Median = ', [num2str(median(IT_EC),' %.2f') ' +/- ' num2str(stdErr, ' %.2f')]};
ylabel('Number of Cells');
xlabel('Correlation');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
yl = ylim;
text(-0.5, yl(2)*0.9, str, 'Color', 'k', 'FontSize', 14, 'FontWeight', 'bold');
filename = [diskPath filesep task filesep 'EncCorr_ITAx&Reacneurons'];
% print(f, filename, '-dpng', '-r0')


%% 

t1 = find(p_cc < 0.05);
t1s = find(p_cc_s < 0.05);
tt = unique([t1 t1s]); % sig corr to scrn or reac

% find the cases where they share and pick the lower value

sigRD = reacData(tt, :);

f = figure; clf;
hold on
IT_EC = sigRD.encCorr(sigRD.reactivated & sigRD.axisTuned);
IT_SC = sigRD.scrnCorr(sigRD.reactivated & sigRD.axisTuned);

% histogram(IT_SC, 8, 'FaceColor', [0.4940 0.1840 0.5560]); stdErr = std(IT_SC)/sqrt(length(IT_SC));
% filename = [diskPath filesep task filesep 'ScrnCorr_ITAx&Reacneurons_SigCorr'];
% str = {'Median = ', [num2str(median(IT_SC),' %.2f') ' +/- ' num2str(stdErr, ' %.2f')]};

histogram(IT_EC, 6, 'FaceColor', [0.4940 0.1840 0.5560]); stdErr = std(IT_EC)/sqrt(length(IT_EC));
filename = [diskPath filesep task filesep 'EncCorr_ITAx&Reacneurons_SigCorr'];
str = {'Median = ', [num2str(median(IT_EC),' %.2f') ' +/- ' num2str(stdErr, ' %.2f')]};

ylabel('Number of Cells');
xlabel('Correlation');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
yl = ylim;
text(0.3, yl(2)*0.9, str, 'Color', 'k', 'FontSize', 14, 'FontWeight', 'bold');
print(f, filename, '-dpng', '-r0')



%% Scatter plot

% enc_r = cell2mat(cellfun(@(x) median(zscore(x)), full_EncResponses(:, 1), 'UniformOutput', false));
% im_r = cell2mat(cellfun(@(x) median(zscore(x)), full_CRResponses(:, 1), 'UniformOutput', false));
% 
% % enc_r = cell2mat(cellfun(@(x) median(x), full_EncResponses(:, 1), 'UniformOutput', false));
% % im_r = cell2mat(cellfun(@(x) median(x), full_CRResponses(:, 1), 'UniformOutput', false));
% 
% reac_ids = find(reacData.reactivated(:));
% nonreac_ids = find(~reacData.reactivated(:));
% 
% f= figure; 
% hold on
% scatter(enc_r(~reacData.reactivated & reacData.axisTuned), im_r(~reacData.reactivated & reacData.axisTuned))
% scatter(enc_r(reacData.reactivated & reacData.axisTuned), im_r(reacData.reactivated & reacData.axisTuned), '*')
% 
% xlim([-0.5 0.5])
% ylim([-0.5 0.5])
% plot([-0.5 0.5], [-0.5 0.5])



%% cdfs



e_cdf = 0;
f = figure;
hold on

c1 = cc_all(~ovrlap);
c2 = cc_all(ovrlap);

% c1 = cc_all(~ovrlap | ~ovrlap_ax);
% c2 = cc_all(ovrlap & ovrlap_ax);


% c1 = cc_all(~ovrlap & ovrlap_ax);
% c2 = cc_all(ovrlap & ovrlap_ax);


if e_cdf
    % % ecdf with bounds - try to make this work
    e1 = ecdf(c1, 'Bounds', 'on'); % ortho
    e2 = ecdf(c2, 'Bounds', 'on');% pref
    grid on
%     plot(e1, x1, 'LineWidth', 2);
%     plot(e2, x2, 'LineWidth', 2);
    
else
    cd1 = cdfplot(c1);
    cd2 = cdfplot(c2);
    cd1.LineWidth = 2;
    cd2.LineWidth = 2;

    [c_p, x_p, ~, ~, ~] = cdfcalc(c1);
    [c_o, x_o, ~, ~, ~] = cdfcalc(c2);
end

[h, p] = kstest2(c1, c2);


% ylim([0 15])
y_lim = ylim;

if e_cdf
    lgnd = legend({'Ortho axis','', '', 'Preferred axis','', ''});
    filename = [diskPath filesep task filesep 'ECDFFRCorr_EncvsCR'];

else
    if ITOnly
        lbl = 'ITOnly';
    else
        lbl = 'AllCells';
    end
    
    lgnd = legend({'Non-Reactivated','Reactivated'});
    filename = [diskPath filesep task filesep 'CDFFRCorr_EncvsCR_' lbl];

end

% if exist('ovrlap', 'var')
%     filename = [filename '_reactivatedCells'];
% end


xlabel('x = Correlation value');
% ylabel('No of neurons');

if strcmp(task, 'Object_Screening')
    xlim([-0.2 1]);
    x_lim = xlim;

    text(0.15*0.80, y_lim(2)*0.88, ['p = ' num2str(p)], 'FontSize', 14,'FontWeight', 'bold');
lgnd.Position = [0.327678571428572,0.646428571428573,0.224107142857143,0.08452380952381];
    title({'Correlation of projection value vs firing rate', 'Preferred and Orthogonal axes', 'Viewing'})
    
elseif strcmp(task, 'Recall_Task')
    x_lim = xlim;
lgnd.Position = [0.170982145837375,0.550000001490118,0.291071422610964,0.110714282734053];

if ITOnly
    title({'Correlation of Responses', 'Viewing and Imagination', 'IT Cells'})
else
    title({'Correlation of Responses', 'Viewing and Imagination', 'All Cells'})

end
    text(x_lim(1)*0.80, y_lim(2)*0.88, ['p = ' num2str(p)], 'FontSize', 14,'FontWeight', 'bold');

end

set(gca, 'FontSize', 14, 'FontWeight', 'bold');
% print(f, filename, '-dpng', '-r0')


%% comparing manual and auto reac

% % which cells did both
% check_ovrlap = ismember(reac_names(:, 1), RC_nm(:, 1));
% % make sure channel numbers are the same too
% for i = 1:length(check_ovrlap)
%     
%     if check_ovrlap(i) == 1
%         rel_idx = find(reac_names(i, 1) == RC_nm(:, 1));
%         
%         if ~isequal(RC_nm(rel_idx, 2), reac_names(i, 2))
%             check_ovrlap(i) = false;
%         end
%     end
%     
% end
% 
% chosenCells = im_active_cells(check_ovrlap, :);
% nonchosenCells = im_active_cells(~check_ovrlap, :);
% % responses = responses(ovrlap, :);
% % psths = psths(ovrlap, :);




