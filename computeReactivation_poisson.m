% Script to compute reactivation via poisson burst metric
% This will allow for assessment of reactivation in individual trials.
% create a folder of specific react cells I want it to capture
setDiskPaths

taskCodePath = [boxPath filesep 'RecallTaskVarun'];
addpath(taskCodePath)
setTTLCodes

% set session list
task = 'Recall_Task';
sessID = Utilities.sessionListAllTasks(task, false);
sessID{2} = ['Recall_Task' filesep 'P76CS' filesep 'Recall_Session_2_20210925']; % want the rasters so using recall session
sessID = sessID';

includeScreening = true;
ITOnly = true;
incAnova = false;

im_active_cells = {};
full_strctCELL = {};
full_ITResp = {};
computeCorr = true;

numReacTrialsPerStim = 3;
numReacTrialsPerCell = 1;

reacTrials = {};

manual_curation = false;
trials_Reac_AllCells = [];

for ss = 1:length(sessID)
    
    tic
%     if includeScreening && ss == 2 % this session didn't have screening
%         continue
%     end
    load([sessID{ss} filesep 'RecallData_NoFreeRec.mat'])
    timelimits_CR = RecallData.offsetTones; 

    % should potentially load responses anyway - for correlation
    % computation
%     if ss ~= 2
%         load([sessID{ss} filesep 'ITResponses.mat']) % strctResp
%     end
    if includeScreening && ss ~=2
        load([sessID{ss} filesep 'PsthandResponses.mat']);
        load([sessID{ss} filesep 'ITResponses.mat']) % strctResp
        full_ITResp = cat(1, full_ITResp, struct2cell(strctResp')');
    end
    
    strctCELL = struct2cell(strctCells')';  
    CRTimeCourse = RecallData.CRTimeCourse;
    EncTimeCourse = RecallData.EncodingTimeCourse;
    
    if ITOnly
        IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
        %         full_strctCELL = strctCELL;
        strctCELL = strctCELL(IT_Cells, :);
        strctCells = strctCells(IT_Cells);
        
        %         BTimeCourse = BTimeCourse(IT_Cells, :);
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
    
    num_cells(ss) = length(strctCells);

    
    for cellIndex = l(strctCells)
        
        addCell = false;
        exCell = CRTimeCourse{cellIndex, 1};
        
%         exCell = exCell(RecallData.CROrder == 2, :); % for testing only
        trials_Reac = nan(size(exCell, 1), 1);
        
        % use reactivation per trial
        for it = 1:size(exCell, 1)
            
            times = find(exCell(it, :) == 1); % should use cutoffs here probably
           
            
            startT = timelimits_CR(1);
            endT = timelimits_CR(1)+5000;
            
            times = find(exCell(it, startT:endT) == 1); % should use cutoffs here probably
            avgSpikRate = sum(times >= startT & times <= endT)/(endT-startT);
            Anchor = 300; % increases rate parameters for poisson
            MinSpInBurst = 5;
            Signif = 0.01;
            MaxXT = 250;
            [b, e, s] = Utilities.p_burst_varun(times, startT, endT, Anchor, MinSpInBurst, Signif, avgSpikRate, MaxXT);
            
            if ~isempty(b)
                %         train = times(times > -timelimits(1)*1e3);
                train = times;%(times > -timelimits(1)*1e3);
                
                % in cases with multiple burtst take one with max surprise
                if length(s) > 1
                    [~, pos] = max(s);
                    stamps(it, 1) = train(b(pos));
                    stamps(it, 2) = train(e(pos));
                    stamps(it, 3) = s(pos);
                else % else take the first
                    stamps(it, 1) = train(b);
                    stamps(it, 2) = train(e);
                    stamps(it, 3) = s;
                end
                trials_Reac(it, 1) = 1;
            else
                trials_Reac(it, 1) = 0;
            end
   
        end
        numRespTrials = nansum(trials_Reac);
        trials_Reac_AllCells(end+1, :) = [strctCells(cellIndex).Name sum(trials_Reac)];
        
        % validating the burst metric - sort and save the psths along
        % with list of reactivated trials
        reacTrials(end+1, :) = [{strctCells(cellIndex).Name} {trials_Reac} {sessID(ss)}];
        
        % ids of stims where reactivated trials 
        stim_Reac = RecallData.CROrder(logical(trials_Reac));  
        
        % if there are n reactive trials for any single stimulus
        cnt = [];
        if numRespTrials >= numReacTrialsPerStim          
            elems = unique(stim_Reac);
            for ii = 1:length(elems)
                cnt(ii, 1) = sum(stim_Reac == elems(ii));
            end
            if sum(cnt >= numReacTrialsPerStim) ~= 0
%                 im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
                addCell = true;
            end        
        end
        
        if incAnova
            cr_psth = CRTimeCourse{cellIndex, 1}(:, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000);
            cr_trials = mean(cr_psth, 2);
            
            group = RecallData.CROrder;
            p_anov = anova1(cr_trials, group, 'off');
        else 
            p_anov = nan;
        end
        
        if addCell || p_anov < 0.05
             im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
        end
        
%         if numRespTrials > numReacTrialsPerCell %(size(exCell, 1)/length(unique(RecallData.CROrder)))/2 % very liberal criteria
%             
%             if incAnova
%                 cr_psth = CRTimeCourse{cellIndex, 1}(:, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000);
%                 cr_trials = mean(cr_psth, 2);
%                 
%                 group = RecallData.CROrder;
%                 p_anov = anova1(cr_trials, group, 'off');
%                 if p_anov < 0.05
%                     im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
%                 end
%             else
%                 im_active_cells = cat(1, im_active_cells, strctCELL(cellIndex, :));
%                 
%             end
%             
%             
%         end
        
        if computeCorr
            % for correlation computation
            if ss == 1
                idx = cellIndex;
            else
                idx = sum(num_cells(1:ss-1)) + cellIndex;
            end
            spikeWinEnc = 500;
%             spikeWinEnc = 267;
            rLUsed = false;
            % grab the encoding responses - like in plotPerCellRasterandbarPerStim
            if exist('responses', 'var') && ~isnan(responses{cellIndex, 2}) && ss ~=2
                start_e = RecallData.offsetEnc(1) + floor(responses{cellIndex, 2});
                till_e = start_e + spikeWinEnc;
                
                start_s = floor(responses{cellIndex, 2}) + 170;
                till_s = start_s + 267;
                
                 rLUsed = true;
            else
                start_e = RecallData.offsetEnc(1)+100;
                till_e = start_e + spikeWinEnc;
                
                start_s = 170+100;
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
            
        end
        
        
    end
    
    full_strctCELL = cat(1, full_strctCELL, strctCELL);
    fprintf('Finished for session %d', ss)
    fprintf('\n')
    toc
    
    
end
full_strctCells = cell2struct(full_strctCELL, fieldnames(strctCells), 2);

% names = cell2mat(im_active_cells(:, 1));
% find(names == 1609)

fprintf('Anchor: %d \n', Anchor)
fprintf('Signif: %d \n', Signif)
fprintf('MinSpikesPerBurst: %d \n', MinSpInBurst)
fprintf('numReactrialsPerStim: %d \n', numReacTrialsPerStim)
fprintf('Window for Enc spikes: %d \n', spikeWinEnc)
fprintf('Response Latency used: %d \n', rLUsed)
fprintf('Screening included?: %d \n', includeScreening)
fprintf('Anova after Poisson?: %d \n', incAnova)
fprintf('Number of Reac Cells: %d \n', size(im_active_cells, 1))



%% compute significance of correlation - per cell
n_reps = 100;
tic
for cellIndex = 1:length(full_strctCells)
    
    enc_s = full_EncResponses{cellIndex, 1};
    CR_s = full_CRResponses{cellIndex, 1};
    
    
    % compute correlation
    cc_boot = nan(n_reps, 1);
    for itr = 1:n_reps
        cc_boot(itr, 1) = corr(Utilities.Shuffle(enc_s), CR_s);
%                 cc_boot(itr, 1) = corr(Utilities.Shuffle(enc_s), CR_s, 'Type', 'Spearman');
    end
    
    cc_all(cellIndex) = corr(enc_s, CR_s);
%         cc_all(cellIndex) = corr(enc_s, CR_s, 'Type', 'Spearman');
    p_cc(cellIndex) = sum(cc_boot > cc_all(cellIndex))/length(cc_boot);
    
    if includeScreening && ~isempty(full_ScrnResponses{cellIndex, 1})
        scrn_s = full_ScrnResponses{cellIndex, 1};
        cc_boot_s = nan(n_reps, 1);
        for itr = 1:n_reps
%                         cc_boot_s(itr, 1) = corr(Utilities.Shuffle(scrn_s), CR_s, 'Type', 'Spearman');
            cc_boot_s(itr, 1) = corr(Utilities.Shuffle(scrn_s), CR_s);
        end
        cc_all_s(cellIndex) = corr(scrn_s, CR_s);
%                 cc_all_s(cellIndex) = corr(scrn_s, CR_s, 'Type', 'Spearman');
        p_cc_s(cellIndex) = sum(cc_boot_s > cc_all_s(cellIndex))/length(cc_boot_s);
    end
    
    
end

fn = fieldnames(strctCells);
RC = cell2struct(im_active_cells, fn(1:6), 2);
toc


%% pulling out various combinations of cells - making reacData
tic
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

fname = ['ReactivationCorrelationData_Poisson_' num2str(Signif)];
if includeScreening
    reacData = table(full_names(:, 1), full_names(:, 2), full_strctCELL(:, 4), cc_all', p_cc', cc_all_s', p_cc_s', ovrlap, ovrlap_resp, ovrlap_ax);
    reacData.Properties.VariableNames = {'Name', 'Channel', 'Area', 'encCorr', 'encCorrSig', 'scrnCorr', 'scrnCorrSig', 'reactivated', 'responsive', 'axisTuned'};
    fname = [fname '_wScrn'];
else
    reacData = table(full_names(:, 1), full_names(:, 2), full_strctCELL(:, 4), cc_all', p_cc', ovrlap, ovrlap_resp, ovrlap_ax);
    reacData.Properties.VariableNames = {'Name', 'Channel', 'Area','encCorr', 'encCorrSig', 'reactivated', 'responsive', 'axisTuned'};
end
if ITOnly
    fname = [fname '_ITOnly'];
end

save([diskPath filesep task filesep 'ReacDataPoisson_' num2str(numReacTrialsPerStim) 'ReacTrPerStim_MinSPB' num2str(MinSpInBurst)], 'reacData', 'reacTrials') 

toc

%% load in axis tuned neurons
tic
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


toc







