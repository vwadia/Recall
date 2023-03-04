%% Script to compute spike field coherence
%
% Try using Juri's code (fieldtrip)
% Try using Ueli's SFCPackage
% Try to write it in house with similar steps
%
% But most important: Learn to prep the data correctly for my sessions - do
% you get coherence spuriously for Im because there's no stim?
%
% LFP prep:
%     Per Session:
%         -- handled in defineusableclusters
%         -- also handled in Ueli's getLFPofTrial function (can read
%         directly from CSC so faster than converting them all to mat
%         files)
%         filter out line noise
%         downsample
%         --
%         high pass filter 1hz (can use Ueli's filterLFPoftrial function)
%
%
% Spike prep:
%     Per Session:
%         Make rasters at 1000hz and separate trials
%   Then can use Jonathan's function for wavelet use or Ueli's for
%   multitaper stuff and compare


% % Ueli's flow
% GetSFC_demo - main func
%     Does all spike prep (gets trialON/OFF times and spike times)
%     precomputeLFPpower - lfp prep
%         getLFPofTrial - grabs each trials lfp (from .ncs file directly), removes line noise and downsamples to 1000Hz
%         filterLFPofTrial - bandpass filter
%         (Now data for each channel is stored as 1xn_trials cell array with lfp of each trial in a cell)
%
%     staExtract_prepare - Prepares data for STA comp (further downsampling)
%     staExtract - gets STA along with it's power spectrum via multitaper method
%     spikeTriggerdSpectrum - ??
%
%     Uses traces for STA and output of sTS to compute SFC

% my pipeline
%     Use Ueli's flow to prep spike data
%     and lfp data (getLFPofTrial)
%
%     Feed into Jonathan's ppc function

% compare multitaper/hilbert/morlet phase extraction methods


%% set paths and load data

setDiskPaths
taskCodePath = [boxPath filesep 'recallTaskVarun']; % more events described here
addpath(taskCodePath); setTTLCodes;

addpath(['Code' filesep 'SFCpackage' filesep 'helpers']);
addpath(genpath('osortTextUI'));

task = 'Recall_Task';
rawPath = 'raw';
sortPath = 'sort';
finalPath = 'final';

% list of sessions
if strcmp(task, 'Recall_Task')
   
    sessID = {['Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917'],...
        ['Recall_Task' filesep 'P76CS' filesep 'Recall_Session_2_20210925'],... % use screening session
        ['Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_3_20210927'],...
        ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_1_20220330'],...
        ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_2_20220403'],...
        ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_3_20220405'],...
        ['Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_1_20220728'],...
        ['Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_2_20220731']};
    
end

invalid_Inds = cell(1, 2); % indices of valid cells
cutOffSpikeVal = 50; % need at least 50 spikes per condition


direction = 'ITCell_HIPPlfp';
% direction = 'ITlfp_HIPPCell';

for cond = 1:2
    
    switch cond
        case 1
            condition = 'Screening';
%             sessID = sessID_1;
        case 2
            condition = 'Imagination';
%             sessID = sessID_2;            
    end
    
    for sess = 1:length(sessID) % start small then do all sessions
        
        basePath = [diskPath filesep sessID{sess}];
        
        % load all events and create trials     
        events = getRawTTLs([basePath filesep rawPath filesep 'Events.nev'], 1);
        expON = find(events(:, 2) == EXPERIMENT_ON);
        expOFF = find(events(:, 2) == EXPERIMENT_OFF);
        
        % grab correct screning session
        if sess == 2
            events_supp = getRawTTLs([diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'RecallScreening_Session_2_20210925' filesep rawPath filesep 'Events.nev'], 1);
            expON_supp = find(events_supp(:, 2) == EXPERIMENT_ON);
            expOFF_supp = find(events_supp(:, 2) == EXPERIMENT_OFF);
            % set path to appropriate screening session
            if cond == 1
                basePath = [diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'RecallScreening_Session_2_20210925'];
            end
        end
        
        switch sess
            case 1
                events1 = events(expON(1):expOFF(1), :); % 1892 trials instead of 2000 makes data 
                events2 = events(find(events(:, 2) == TRAINING_RECALL_END):expOFF(end), :); % recall task
            case 2 % pain in the ass session
                events1 = events_supp(expON_supp(1):expOFF_supp(1), :); % screening
                events2 = events(find(events(:, 2) == TRAINING_RECALL_END):expOFF(end), :); % recall task
            case 3
                events1 = events(expON(end):expOFF(end), :); % screening
                events2 = events(find(events(:, 2) == TRAINING_RECALL_END):expOFF(1), :); % recall task
            case 4
                events1 = events(expON(end):expOFF(end), :); % screening
                trainEnd = find(events(:, 2) == TRAINING_RECALL_END);
                events2 = events(trainEnd(end):expOFF(1), :); % recall task
            case 5
                events1 = events(expON(end):expOFF(end), :); % screening
                events2 = events(find(events(:, 2) == TRAINING_RECALL_END):expOFF(1), :); % recall task
            case 6
                events1 = events(expON(end):expOFF(end), :); % screening
                trainEnd = find(events(:, 2) == TRAINING_RECALL_END);
                events2 = events(trainEnd(end):expOFF(2), :); % recall task
            case 7
                events1 = events(expON(end):expOFF(end), :); % screening
                events2 = events(find(events(:, 2) == TRAINING_RECALL_END):expOFF(1), :); % recall task
            case 8
                events1 = events(expON(end):expOFF(end), :); % screening
                events2 = events(find(events(:, 2) == TRAINING_RECALL_END):expOFF(1), :); % recall task                
        end
        
        events1MS = events1*1e3;
        events2MS = events2*1e3;
        
        % extract all neurons
        [strctCells, dupCells] = Utilities.extractCells([basePath filesep sortPath filesep finalPath], basePath);
        strctCELL = struct2cell(strctCells');
        strctCELL = strctCELL';
        
        % narrow down analysis regions
        if strcmp(direction, 'ITCell_HIPPlfp')
            %         IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
            IT_Cells = cellfun(@(x) strcmp(x, 'LFFA'), strctCELL(:, 4));
            strctCells = strctCells(IT_Cells);
            
            % define lfp channels
            % LH = 25:32; RH = 57:64; LIT = 209:216; RIT = 217:224;
            LChans = 25:32;
            RChans = 57:64;
        elseif strcmp(direction, 'ITlfp_HIPPCell')
            Hipp_Cells = cellfun(@(x) strcmp(x, 'RH') || strcmp(x, 'LH'), strctCELL(:, 4));
            strctCells = strctCells(Hipp_Cells);
            
            % define lfp channels
            % LH = 25:32; RH = 57:64; LIT = 209:216; RIT = 217:224;
            LChans = 209:216;
            RChans = 217:224;
        end
        
        % define time periods
        % screening - (stim on - 100ms) so you don't miss 1st trial
        % note that resolution is in microseconds so 100ms = 1e5 us
        periods1 = [events1(events1(:, 2) == IMAGE_ON, 1) - 1e5 events1(events1(:, 2) == IMAGE_ON, 1)];
        
        % recall
        periods2 = [events2(events2(:, 2) == TONE_1, 1) events2(events2(:, 2) == TONE_2, 1);...
            events2(find(events2(:, 2) == TONE_2), 1) events2(find(events2(:, 2) == TONE_2)+1, 1)];
        
        if strcmp(condition, 'Imagination')
            timeFrom = periods2(:, 1);
            timeTo = periods2(:, 2);
            bin_size = 200;
            Ord = 1:length(timeFrom);
            use_both_offsets = 1;
            
            offset = [0 0];
            
        elseif strcmp(condition, 'Screening')
            timeFrom = periods1(:, 1);
            timeTo = periods1(:, 2);
            
            bin_size = mean(timeTo-timeFrom)*1e-4; % because timeFrom and To are in us
            Ord = 1:length(timeFrom);
            use_both_offsets = 1;
            
            offset = [170 530]; % in ms - should this be us??
        end
        
        data_spike = [];
        % make rasters & collect spike data
        for cellIndex = l(strctCells)
            [psths{cellIndex, 1}, psths{cellIndex, 2},  psths{cellIndex, 3}]...
                = Utilities.makeRastersFromStrctCells(strctCells(cellIndex), offset, use_both_offsets, Ord, bin_size, timeFrom*1e-3, timeTo*1e-3);  % note conversion to ms
            
            sDat = psths{cellIndex, 1}'; % samples x trials
            
            sDat = reshape(sDat, [size(sDat, 1) 1 size(sDat, 2)]);
            
            if ~isempty(data_spike)
                if size(data_spike, 1) > size(sDat, 1)
                    sDat = padarray(sDat, size(data_spike, 1) - size(sDat, 1), 'post');
                elseif size(data_spike, 1) < size(sDat, 1)
                    sDat = sDat(1:size(data_spike, 1), :, :);
                end
            end
            
            % put into correct shape for Jonathan's function
            data_spike = cat(2, data_spike, sDat);
        end
        
        
        % ppc parameters
        low_freq = 2;
        high_freq = 100;
        n_freq = high_freq - low_freq;
        
        
        % parameters for getting lfp of trials
        useNotchFilter=0;
        replaceSpikes=0;
        if sess == 2 && strcmp(condition, 'Screening')
            basepathLFP=[diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'RecallScreening_Session_2_20210925' filesep rawPath];
        else
            basepathLFP=[basePath filesep rawPath];
        end
        replaceMethod=1; %1=spline, 2=cubic, 3=linear, 4=mean spike subtraction
        cancelSpikePhaseShift=0;
        FsDownReq=1000; % requested sampling rate.
        data_lfp = [];
        
        % grab lfp data
        for channel = LChans%[LChans RChans]
            dfLP = [];
            
            % data is returned as a 1 x ntrials cell array
            [dataRaw, dataRawDownsampled, isContinuous,periodsExtracted, validTrials, FsDown,Fs, maxRange] = getLFPofTrial(basepathLFP, channel, timeFrom, timeTo, useNotchFilter, FsDownReq );
            
            %filter
            filterMode=9; %LFP
            [dataFilteredLowPass, filterSettings1] = filterLFPofTrial( dataRawDownsampled, filterMode, FsDown );
            
            
            maxLength = cellfun(@(x) length(x), dataFilteredLowPass);
            maxLength = max(maxLength);
            
            % make sure size matches spike data
            if size(data_spike, 1) ~= maxLength
                maxLength = size(data_spike, 1);
            end
            
            for tr = 1:length(dataFilteredLowPass)
                if length(dataFilteredLowPass{tr}) > maxLength
                    dataFilteredLowPass{tr} = dataFilteredLowPass{tr}(1:maxLength, :);
                elseif length(dataFilteredLowPass{tr}) < maxLength
                    dataFilteredLowPass{tr} = padarray(dataFilteredLowPass{tr}, maxLength - length(dataFilteredLowPass{tr}), 'post');
                end
            end
            
            % put it into correct shape for jonathan's function
            dfLP = cell2mat(dataFilteredLowPass);
            dataRawReshaped = reshape(dfLP, [1 size(dfLP, 1)*size(dfLP, 2)]);
            
            if std(dataRawReshaped) > 0.2 % make sure it was not used as a reference channel
                dfLP = reshape(dfLP, [size(dfLP, 1) 1 size(dfLP, 2)]);
                data_lfp = cat(2, data_lfp, dfLP);
            end
            
        end
        
        % compute ppc for that session and condition
        if strcmp(condition, 'Screening')
            [ppc_Sc{sess}, f_Sc{sess}] = Recall.compute_ppc(data_lfp, data_spike, low_freq, high_freq, n_freq, 'linear', FsDown);
        elseif strcmp(condition, 'Imagination')
            [ppc_Im{sess}, f_Im{sess}] = Recall.compute_ppc(data_lfp, data_spike, low_freq, high_freq, n_freq, 'linear', FsDown);
        end
    end
end

%% plot both conditions separately showing each session

figSc = figure;
hold on
for s = 1:length(sessID)
    
    plot(f_Sc{s}, squeeze(nanmean(nanmean(ppc_Sc{s}))), 'LineWidth', 2);
%     keyboard
end

figIm = figure;
hold on
for s = 1:length(sessID)
    
    plot(f_Im{s}, squeeze(nanmean(nanmean(ppc_Im{s}))), 'LineWidth', 2);
%     keyboard
end


%% all together summary figure
avPPC_Sc = [];
avPPC_Im = [];
for s = 1:length(sessID)
    
    avPPC_Sc = [avPPC_Sc squeeze(nanmean(nanmean(ppc_Sc{s})))];
    avPPC_Im = [avPPC_Im squeeze(nanmean(nanmean(ppc_Im{s})))];
    
end

avPPC_Sc = nanmean(avPPC_Sc, 2);
avPPC_Im = nanmean(avPPC_Im, 2);

sumFig = figure;
hold on
plot(f_Sc{1}, avPPC_Sc, 'LineWidth', 2)
plot(f_Im{1}, avPPC_Im, 'LineWidth', 2)
title('PPC LIT cells - LH lfp')
legend('BottomUp', 'TopDown');
set(gca, 'FontSize', 14, 'FontWeight', 'bold')





























