% Better (cleaner) version of compute_SFC
% Improvements:
%     Gather all data (spike/Lfp) for both conditions *first*
%         Instead of re-extracting cells load them in per session (dummy)
%     Choose only ramp tuned neurons? Or visually responsive neurons?
%     Then based on condition/side/direction - maybe write different script for Hippcell - ITlfp


%% set paths and define sessions
setDiskPaths
taskCodePath = [boxPath filesep 'recallTaskVarun']; % more events described here
addpath(taskCodePath); setTTLCodes;

addpath([diskPath filesep 'Code' filesep 'SFCpackage' filesep 'helpers']);
addpath([diskPath filesep 'Code' filesep 'SFCpackage' filesep 'SFCieeg']);
addpath(genpath('osortTextUI'));

task = 'Recall_Task';
rawPath = 'raw';
sortPath = 'sort';
finalPath = 'final';

% session details - note these are for RECALL (relevant for P76 sess 2)
dirID = {['Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917'];...
    ['Recall_Task' filesep 'P76CS' filesep 'Recall_Session_2_20210925'];... % use screening session
    ['Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_3_20210927'];...
    ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_1_20220330'];...
    ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_2_20220403'];...
    ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_3_20220405'];...
    ['Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_1_20220728'];...
    ['Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_2_20220731']};
sessionID = {'P76CSRec_ReScreen';  'P76CS_RecScreen3'; 'P76CSRec_ReScreen_3'; 'P79CS_ReScreen_1'; 'P79CS_ReScreen_3'; 'P79CS_ReScreen_4'; 'P80CS_ReScreenRecall'; 'P80CS_ReScreecRecall_2'};

Side = 'RFFA';
% Side = 'LFFA';

% need at least 50 spikes per condition  for IT
% change this for other Areas
cutOffSpikeVal = 50; 

% params.run_boot = 'false';
params.run_boot = 'true';

C1_dirID(:, 1) = {[1 5 7 8 14 17:24 25 26 28 34:36 38 40 44:46 50 52 54:56 57 60:64 193:195 197:200 211 213 ],... %P76 RsR 1
    [4 6 7 8 14 21 22 24 41 45 50 64 198 211 213],... P76 Screen 2
    [6:8 14 21 22 28 34 45 50 53 56 59 64 197 200 211 213],... P76 RsR 3
    [1:3 5:10 12:16 18:22 24:28 30:33 36 38 46 49:56 59:62 64 193:200 202 208],... P79 RsR 1, all IT channels seem to have noise (ugh)
    [1:3 5:10 13:22 24 27 28 29 30 35 38 46 49:52 54:56 59 60 193 194 196:199 209:212 214 215 217:219 222],... P79 RsR 2
    [1 3 5:10 14:17 19:22 24 27 28 30 32 35 46 49:52 54 55 59 60 62 193:199 202 205 209:212 214 215 217:222],... P79 RsR 3
    [1 3 9 16:25 31 32 35 37:41 43 45:56 59 206 209:211 214 217 222:224],... P80 RsR 1
    [1:5 7 9:11 16:26 30:32 35 37:39 41:56 58 59 206 209 211 219 222:224]}; % P80 RsR 2

C2_dirID(:, 1) = {[1 5 7 8 14 17:24 25 26 28 34:36 38 40 44:46 50 52 54:56 57 60:64 193:195 197:200 211 213 ],... %P76 RsR 1
    [4 6:8 14 21 22 24 28 34 38 40 41 44 45 50 53 193 195 200 205 211 213],... P76 Recall 2
    [6:8 14 21 22 28 34 45 50 53 56 59 64 197 200 211 213],... P76 RsR 3
    [1:3 5:10 12:16 18:22 24:28 30:33 36 38 46 49:56 59:62 64 193:200 202 208],... P79 RsR 1, all IT channels seem to have noise (ugh)
    [1:3 5:10 13:22 24 27 28 29 30 35 38 46 49:52 54:56 59 60 193 194 196:199 209:212 214 215 217:219 222],... P79 RsR 2
    [1 3 5:10 14:17 19:22 24 27 28 30 32 35 46 49:52 54 55 59 60 62 193:199 202 205 209:212 214 215 217:222],... P79 RsR 3
    [1 3 9 16:25 31 32 35 37:41 43 45:56 59 206 209:211 214 217 222:224],... P80 RsR 1
    [1:5 7 9:11 16:26 30:32 35 37:39 41:56 58 59 206 209 211 219 222:224]}; % P80 RsR 2

%% load in spike data for both conditions and find valid neurons

% load in cells 
%  load([diskPath filesep task filesep 'AllITCells_500stim_Im.mat']);
%  load([diskPath filesep task filesep 'AllRespITCells_500stim_Im.mat']);
load([diskPath filesep task filesep 'AllITCells_500Stim_Im_SigRamp.mat']);
strctCELL = struct2cell(strctCells');
strctCELL = strctCELL';
 
tic
% For a given session 
for sess = 1:length(dirID)
    
    % this fucking second session - loading in all cells won't work for
    % screening
%     if sess == 2
%         basePath = [diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'RecallScreening_Session_2_20210925'];
%     else
    basePath = [diskPath filesep dirID{sess}];
%     end
    % carve out relevant neurons from screening
    relevantCells = cellfun(@(x,y) strcmp(x, Side) && strcmp(y, sessionID{sess}), strctCELL(:, 4), strctCELL(:, 8), 'UniformOutput', false);
    relevantCells = cell2mat(relevantCells);
    sess_strctCells = strctCells(relevantCells);
    sess_psths = psths(relevantCells, :);
    sess_responses = responses(relevantCells, :);
    timelimits_scrn = [-0.17 0.53];
    timelimits_im = [-1 6.5];
    offset_scrn = 100;
    offset_im = 5000;
    
    sData = cell(length(sess_strctCells), 1);
    % grab relevant spike data for screening
    for cellIndex = l(sess_strctCells)
        sData{cellIndex, 1} = sess_psths{cellIndex, 1}(:, -timelimits_scrn(1)*1e3-offset_scrn:-timelimits_scrn(1)*1e3)';
    end
    
    % load in Im data
    ImStrct = load([basePath filesep 'RecallData_NoFReeRec']);
    ImStrct.strctCELL = struct2cell(ImStrct.strctCells');
    ImStrct.strctCELL = ImStrct.strctCELL';
    relevantCells = cellfun(@(x) strcmp(x, Side), ImStrct.strctCELL(:, 4)); % cells in the correct area
    
    ImStrct.strctCells = ImStrct.strctCells(relevantCells);
    ImStrct.RecallData.CRTimeCourse = ImStrct.RecallData.CRTimeCourse(relevantCells, :);
    ImStrct.RecallData.EncodingTimeCourse = ImStrct.RecallData.EncodingTimeCourse(relevantCells, :);
    
    % correctly pick out imagination trial data for responsive/ramp tuned neurons
    for cellIndex = l(sess_strctCells)
        cntrl = 1;
        
        % manually matching cells from morning and afternoon for this
        % annoying session
        if sess == 2
            matches = [2 2; 3 4; 4 3];
            if ~isempty(find(matches(:, 1) == cellIndex))
                disp(matches(find(matches(:, 1) == cellIndex), 2))
                sData{cellIndex, 2} = ImStrct.RecallData.CRTimeCourse{matches(find(matches(:, 1) == cellIndex), 2), 1}(:, -timelimits_im(1)*1e3:-timelimits_im(1)*1e3+offset_im)';
            else
                sData{cellIndex, 1} = [];
            end
        else
            for c_idx = l(ImStrct.strctCells)



                if isequal(sess_strctCells(cellIndex).Name, ImStrct.strctCells(c_idx).Name) &&...
                        isequal(sess_strctCells(cellIndex).spikeTimeStamps, ImStrct.strctCells(c_idx).spikeTimeStamps)

                    assert(cntrl == 1, 'Cell Pairing is wrong!');
                    % sData needs to be cellIndex and EncTimeCourse needs to be c_idx
                    sData{cellIndex, 2} = ImStrct.RecallData.CRTimeCourse{c_idx, 1}(:, -timelimits_im(1)*1e3:-timelimits_im(1)*1e3+offset_im)';
                    cntrl = cntrl + 1;
                end
            end
        end
    end
    if sess == 2
        sData = sData(2:end, :);
    end

    % Check that cells spike enough and grab only valid ones
    % Do I need to do this if it's PPC?
    valid_cell = zeros(length(sess_strctCells), 1);
    if sess ~= 2
        assert(isequal(length(valid_cell), size(sData, 1)));
    end
    for cellIndex = 1:(size(sData, 1))
        
        if sum(sum(sData{cellIndex, 1})) > cutOffSpikeVal &&...
                sum(sum(sData{cellIndex, 2})) > cutOffSpikeVal
            
            valid_cell(cellIndex, 1) = 1;
        end
    end
    
    sData = sData(logical(valid_cell), :);
    
    %% event info lfp
    % load all events and create trials
    events = getRawTTLs([basePath filesep rawPath filesep 'Events.nev'], 1);
    expON = find(events(:, 2) == EXPERIMENT_ON);
    expOFF = find(events(:, 2) == EXPERIMENT_OFF);

    % grab correct screning session
    if sess == 2
        events_supp = getRawTTLs([diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'RecallScreening_Session_2_20210925' filesep rawPath filesep 'Events.nev'], 1);
        expON_supp = find(events_supp(:, 2) == EXPERIMENT_ON);
        expOFF_supp = find(events_supp(:, 2) == EXPERIMENT_OFF);
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

    % define time periods
    % screening - (stim on - 100ms) so you don't miss 1st trial
    % note that resolution is in microseconds so 100ms = 1e5 us
    periods1 = [events1(events1(:, 2) == IMAGE_ON, 1) - 1e5 events1(events1(:, 2) == IMAGE_ON, 1)];

    % recall
    periods2 = [events2(events2(:, 2) == TONE_1, 1) events2(events2(:, 2) == TONE_2, 1);...
        events2(find(events2(:, 2) == TONE_2), 1) events2(find(events2(:, 2) == TONE_2)+1, 1)];


    %% compute ppc using Jonathans function for each condition
    for cond = 1:2%size(sData, 2)
         % define lfp channels
        % HIPP
        LChans = 25:32;
        RChans = 57:64;
        % MFC
%         LChans = 9:16;
%         RChans = 41:48;
        % OFC
%         LChans = 193:200;
%         RChans = 201:208; % P76 doesn't have this

        switch cond
            case 1
                condition = 'Screening';
%                 condition = 'Encoding';
                validChans = C1_dirID{sess};
                
            case 2
                condition = 'Imagination';
                validChans = C2_dirID{sess};
        end
        % choose only channels with cells
        LChans = LChans(ismember(LChans, validChans));
        RChans = RChans(ismember(RChans, validChans));
        % put data into correct format for Jonathan's function
        data_spike = [];
        for cellIndex = 1:size(sData, 1)
            
            sD = sData{cellIndex, cond}; 
            sD = reshape(sD, [size(sD, 1) 1 size(sD, 2)]);
            if ~isempty(data_spike)
                if size(data_spike, 1) > size(sD, 1)
                    sD = padarray(sD, size(data_spike, 1) - size(sD, 1), 'post');
                elseif size(data_spike, 1) < size(sD, 1)
                    sD = sD(1:size(data_spike, 1), :, :);
                end
            end
            
            % should be samples x neurons x trials
            data_spike = cat(2, data_spike, sD);
        end
        
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
        
       
        % ppc parameters
        low_freq = 2;
        high_freq = 100;
        n_freq = high_freq - low_freq;
        
        
        % parameters for getting lfp of trials
        useNotchFilter=1;
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
        if strcmp(Side, 'LFFA')
            chans = LChans;
        elseif strcmp(Side, 'RFFA')
            chans = RChans;
        end
        dataRawReshaped = [];
        chan_ind = 1;
        dfLP = cell(1, length(chans));

        % parfor this?
%         for channel = chans
        parfor chan_ind = 1:length(chans)
            channel = chans(chan_ind);
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
            
            % put data into correct shape for jonathan's function
            dfLP{chan_ind} = cell2mat(dataFilteredLowPass);
%             dataRawReshaped = reshape(dfLP, [1 size(dfLP, 1)*size(dfLP, 2)]);
            dataRawReshaped = [dataRawReshaped; reshape( dfLP{chan_ind}, [1 size( dfLP{chan_ind}, 1)*size(dfLP{chan_ind}, 2)])];

%             chan_ind = chan_ind+1;
        end
        
        chan_stddevs = std(dataRawReshaped');
        if mean(chan_stddevs)/min(chan_stddevs) > 2 % means one is reference channel
            [~, ref] = min(chan_stddevs);
             dfLP(ref) = [];
        end
        
        for idx = 1:length(dfLP)
            dfLP_ofChan = reshape(dfLP{idx}, [size(dfLP{idx}, 1) 1 size(dfLP{idx}, 2)]);
            data_lfp = cat(2, data_lfp, dfLP_ofChan);
        end
        
        % compute ppc for that session and condition
        FsDown = FsDownReq;
        if ~isempty(data_lfp)            
            if strcmp(condition, 'Screening')
                [ppc_Sc{sess}, f_Sc{sess}, ppc_Sc_boot{sess}] = Recall.compute_ppc(data_lfp, data_spike, low_freq, high_freq, n_freq, 'log', FsDown, params.run_boot);
                
            elseif strcmp(condition, 'Imagination')
                [ppc_Im{sess}, f_Im{sess}, ppc_Im_boot{sess}] = Recall.compute_ppc(data_lfp, data_spike, low_freq, high_freq, n_freq, 'log', FsDown, params.run_boot);
            end
        end
%         keyboard
    end
    
end
toc
frq = f_Sc{1};
save([diskPath filesep task filesep 'ppc_log' filesep 'ppc_ITCellHippLFPCellChans_sigRamp_wBootLogScale_Right.mat'], 'ppc_Sc', 'ppc_Im', 'frq', 'ppc_Sc_boot', 'ppc_Im_boot')

if strcmp(Side, 'LFFA')
    
    ppc_Sc = ppc_Sc(4:end);
    f_Sc = f_Sc(4:end);
    ppc_Im = ppc_Im(4:end);
    f_Im = f_Im(4:end);
end

%% summary figure

setDiskPaths
task = 'Recall_Task';
% scale = 'ppc_linear';
scale = 'ppc_log';

if strcmp(scale, 'ppc_linear')
    frq = linspace(2, 100, 98);
    load([diskPath filesep task filesep scale filesep 'ppc_RITCellRHippLFP_allCells.mat']); noSess2 = 1; Side = 'RFFA';
%     load([diskPath filesep task filesep scale filesep 'ppc_RITCellRHippLFP_respCells.mat']); noSess2 = 1; singleCellSess8 = 1; Side = 'RFFA';
%     load([diskPath filesep task filesep scale filesep 'ppc_RITCellRHippLFP_SigRamp.mat']); noSess8 = 1; Side = 'RFFA';
   
%     load([diskPath filesep task filesep scale filesep 'ppc_LITCellLHippLFP_allCells.mat']); Side = 'LFFA';
%     load([diskPath filesep task filesep scale filesep 'ppc_LITCellLHippLFP_respCells.mat']); Side = 'LFFA';
%     load([diskPath filesep task filesep scale filesep 'ppc_LITCellLHippLFP_SigRamp.mat']); Side = 'LFFA';
elseif strcmp(scale, 'ppc_log')
%     load([diskPath filesep task filesep scale filesep 'ppc_ITCellHippLFP_allCells_wBootLogScale_Right.mat']); Side = 'RFFA';
%     load([diskPath filesep task filesep scale filesep 'ppc_ITCellHippLFP_SigRamp_wBootLogScale_Right.mat']); Side = 'RFFA';
%     load([diskPath filesep task filesep scale filesep 'ppc_ITCellHippLFPCellChans_sigRamp_wBootLogScale_Right.mat']); Side = 'RFFA';

%     load([diskPath filesep task filesep scale filesep 'ppc_ITCellHippLFP_allCells_wBootLogScale__Left.mat']); Side = 'LFFA';
%     load([diskPath filesep task filesep scale filesep 'ppc_ITCellHippLFP_SigRamp_wBootLogScale_Left.mat']); Side = 'LFFA';
    load([diskPath filesep task filesep scale filesep 'ppc_ITCellHippLFPCellChans_sigRamp_wBootLogScale_Left.mat']); Side = 'LFFA';

    
end

cols = Utilities.distinguishable_colors(2);
avPPC_Sc = [];
avPPC_Im = [];
for s = 1:length(ppc_Sc)
    if ~isempty(ppc_Sc{s}) && ~isempty(ppc_Im{s})
        if size(ppc_Sc{s}, 1) == 1 || size(ppc_Sc{s}, 2) == 1
            avPPC_Sc = [avPPC_Sc squeeze(mean(ppc_Sc{s}))];
        elseif size(ppc_Im{s}, 1) == 1 || size(ppc_Im{s}, 2) == 1
            avPPC_Im = [avPPC_Im squeeze(mean(ppc_Im{s}))];
        else
            avPPC_Sc = [avPPC_Sc squeeze(mean(mean(ppc_Sc{s})))];
            avPPC_Im = [avPPC_Im squeeze(mean(mean(ppc_Im{s})))];
        end                
    end
end

% old way of plotting std error of mean - just use stdshade instead ------
% avPPC_Sc = mean(avPPC_Sc, 2);
% avPPC_Im = mean(avPPC_Im, 2);
% stdErr_Sc = std(avPPC_Sc)/sqrt(length(avPPC_Sc));
% stdErr_Im = std(avPPC_Im)/sqrt(length(avPPC_Im));
% ------------------------------------------------------------------------

sumFig = figure;
% plot(f_Sc{1}, avPPC_Sc, 'LineWidth', 2)
% plot(f_Im{1}, avPPC_Im, 'LineWidth', 2)
if strcmp(scale, 'ppc_linear')
    hold on
    Utilities.stdshade5(avPPC_Sc', 0.1, cols(1, :), frq');
    Utilities.stdshade5(avPPC_Im', 0.1, cols(2, :), frq');
elseif strcmp(scale, 'ppc_log')
    % hold on removes the log axis....
    s1 = semilogx(frq, mean(avPPC_Sc, 2), 'LineWidth', 2);
    s2 = semilogx(frq, mean(avPPC_Im, 2), 'LineWidth', 2);
end
if strcmp(Side, 'RFFA')
    title('PPC RIT cells - RH lfp')
%     ylim([0 10e-4])

elseif strcmp(Side, 'LFFA')
    title('PPC LIT cells - LH lfp')
%     ylim([0 15e-4])

end
filename = [diskPath filesep task filesep scale filesep 'PPC_' Side 'Cell_' Side(1) 'HLFP_AllSessions_allCells_IncSess2'];
% filename = [diskPath filesep task filesep 'PPC_' Side 'Cell_' Side(1) 'HLFP_AllSessions_AllResp'];
legend('BottomUp', 'TopDown');
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
% print(sumFig, filename, '-dpng', '-r0')


%% per session
setDiskPaths
task = 'Recall_Task';

% scale = 'ppc_linear';
scale = 'ppc_log';

if strcmp(scale, 'ppc_linear')
    frq = linspace(2, 100, 98);
%     load([diskPath filesep task filesep scale filesep 'ppc_RITCellRHippLFP_allCells.mat']); noSess2 = 1; Side = 'RFFA';
%     load([diskPath filesep task filesep scale filesep 'ppc_RITCellRHippLFP_respCells.mat']); noSess2 = 1; singleCellSess8 = 1; Side = 'RFFA';
%     load([diskPath filesep task filesep scale filesep 'ppc_RITCellRHippLFP_SigRamp.mat']); noSess8 = 1; Side = 'RFFA';
%     
%     
%     load([diskPath filesep task filesep scale filesep 'ppc_LITCellLHippLFP_allCells.mat']); Side = 'LFFA';
%     load([diskPath filesep task filesep scale filesep 'ppc_LITCellLHippLFP_respCells.mat']); singleCellSess5 = 1; Side = 'LFFA';
%     load([diskPath filesep task filesep scale filesep 'ppc_LITCellLHippLFP_SigRamp.mat']); singleCellSess5 = 1; Side = 'LFFA';
elseif strcmp(scale, 'ppc_log')
%     load([diskPath filesep task filesep scale filesep 'ppc_ITCellHippLFP_allCells_wBootLogScale_Right.mat']); noSess2 = 1; Side = 'RFFA';
%     load([diskPath filesep task filesep scale filesep 'ppc_ITCellHippLFP_SigRamp_wBootLogScale_Right.mat']); noSess2 = 1; Side = 'RFFA';
    load([diskPath filesep task filesep scale filesep 'ppc_ITCellHippLFPCellChans_sigRamp_wBootLogScale_Right.mat']); Side = 'RFFA'; cellChans = 1;
%     
%     load([diskPath filesep task filesep scale filesep 'ppc_ITCellHippLFP_allCells_wBootLogScale__Left.mat']); Side = 'LFFA';
%     load([diskPath filesep task filesep scale filesep 'ppc_ITCellHippLFP_SigRamp_wBootLogScale_Left.mat']); Side = 'LFFA';
%     load([diskPath filesep task filesep scale filesep 'ppc_ITCellHippLFPCellChans_sigRamp_wBootLogScale_Right.mat']); Side = 'LFFA'; cellChans = 1;

end
%%

cols = Utilities.distinguishable_colors(2);
avPPC_Sc = [];
avPPC_Im = [];
for session = 1:length(ppc_Sc)
    if ~isempty(ppc_Sc{session}) && ~isempty(ppc_Im{session})
        if size(ppc_Sc{session}, 1) == 1 || size(ppc_Sc{session}, 2) == 1
            avPPC_Sc = squeeze(ppc_Sc{session});
        elseif size(ppc_Im{session}, 1) == 1 || size(ppc_Im{session}, 2) == 1
            avPPC_Im = squeeze(ppc_Im{session});
        else
            avPPC_Sc = squeeze(mean(ppc_Sc{session}, 2));
            avPPC_Im = squeeze(mean(ppc_Im{session}, 2));
        end
        f = figure; 
        if strcmp(scale, 'ppc_linear')
            hold on

            Utilities.stdshade5(avPPC_Sc, 0.1, cols(1, :), frq');
            Utilities.stdshade5(avPPC_Im, 0.1, cols(2, :), frq');
        elseif strcmp(scale, 'ppc_log')
           c1_ppc = mean(avPPC_Sc);
           stderr_c1 = mean(c1_ppc)/sqrt(length(c1_ppc));
           c2_ppc = mean(avPPC_Im);
           stderr_c2 = mean(c2_ppc)/sqrt(length(c2_ppc));
           semilogx(frq, c1_ppc, frq, c2_ppc, 'LineWidth', 2)
%            patch('Faces', [frq fliplr(frq)], [c1_ppc + stderr_c1 fliplr(c1_ppc - stderr_c1)], 'FaceColor', [0 0 1],'FaceAlpha', 0.1)
        end
        
%         xticks(frq);
%         set(gca, 'xticklabel', {[]}); % clear the old labels
%         xtik = get(gca,'xtick');
%         xtiklabs = {frq};
        yl = ylim;
        ylim([0 yl(2)])
        if exist('cellChans', 'var')
            filename = [diskPath filesep task filesep scale filesep 'ppc_' Side 'Cell_' Side(1) 'HLFP_Sess' num2str(session) '_allCells_cellChans'];
%             filename = [diskPath filesep task filesep scale filesep 'ppc_' Side 'Cell_' Side(1) 'HLFP_Sess' num2str(session) '_sigRamp_cellChans'];            
        else
            filename = [diskPath filesep task filesep scale filesep 'ppc_' Side 'Cell_' Side(1) 'HLFP_Sess' num2str(session) '_allCells'];
%             filename = [diskPath filesep task filesep scale filesep 'ppc_' Side 'Cell_' Side(1) 'HLFP_Sess' num2str(session) '_sigRamp'];
        end
        legend('BottomUp', 'TopDown');
        set(gca, 'FontSize', 14, 'FontWeight', 'bold')
%         print(f, filename, '-dpng', '-r0')
%         close all
    end
    
end
%% Per cell 

chosenSess = 5;

ppc_ScSess = ppc_Sc{chosenSess};
ppc_ImSess = ppc_Im{chosenSess};

for cellIndex = 1:size(ppc_ScSess, 1)
    figure;
    hold on
    cellPPC_Sc = squeeze(mean(ppc_ScSess, 2));%, 'LineWidth', 2);
    cellPPC_Im = squeeze(mean(ppc_ImSess, 2));%, 'LineWidth', 2);
    
    if chosenSess == 8 && strcmp(Side, 'RFFA')
        plot(cellPPC_Sc)
        plot(cellPPC_Im)
    else
        plot(cellPPC_Sc(cellIndex, :)','LineWidth', 2)
        plot(cellPPC_Im(cellIndex, :)','LineWidth', 2)
    end
    yl = ylim;
%     ylim([0 yl(2)])
    keyboard
    
end

%% messing with surrogates

thetaRange = f_Im{1} > 3 & f_Im{1} < 8; % 3 - 8Hz
BetaRange = f_Im{1} > 13 & f_Im{1} < 30; % 13 - 30Hz
lowGammaRange = f_Im{1} > 30 & f_Im{1} < 60; % 30 - 60Hz

for session = 1:length(ppc_Im)
    for neuron = 1:size(ppc_Im{session}, 1)
        for channel = 1:size(ppc_Im{session}, 2)

            ppcIndiv = squeeze(ppc_Im{session}(neuron, channel, thetaRange));
            ppcDistRange = squeeze(ppc_Im_boot{session}(neuron, channel, thetaRange, :));

            for fr = 1:sum(thetaRange)
                ppcDist = ~isnan(ppcDistRange(fr, :));
                ppcDist = ppcDistRange(fr, ppcDist);
                [mu, sig] = normfit(ppcDist);

                % mean subtract
                ppcSub = (ppcIndiv - mu)/sig;

                sigPPC = ppcIndiv(ppcSub > 2);
                sigF = f_Im{1}(thetaRange);
                sigF = sigF(ppcSub > 2);
            end
        end
    end
end

%% Cluster based parametric statistics
setDiskPaths

% load([diskPath filesep 'Recall_Task' filesep 'ppc_ITCellHippLFP_allCells_wBootLogScale_Right.mat']);
% load([diskPath filesep 'Recall_Task' filesep 'ppc_ITCellHippLFP_allCells_wBootLogScale_Left.mat']);




addpath([diskPath filesep 'Code' filesep 'fieldtrip-20200409']) % Jonathan's version
sessStats = {};
for session = 1:length(ppc_Sc)

    if ~isempty(ppc_Sc{session})

        ppc_C1 = ppc_Sc{session};
        n_neurons_C1 = size(ppc_C1, 1);
        n_channels_C1 = size(ppc_C1, 2);
        n_freq_C1 = size(ppc_C1, 3);

        ppc_C1 = reshape(ppc_C1, [n_neurons_C1*n_channels_C1 n_freq_C1]);

        ppc_C2 = ppc_Im{session};
        n_neurons_C2 = size(ppc_C2, 1);
        n_channels_C2 = size(ppc_C2, 2);
        n_freq_C2 = size(ppc_C2, 3);

        ppc_C2 = reshape(ppc_C2, [n_neurons_C2*n_channels_C2 n_freq_C2]);

        frqs = frq;

        sessStats{session} = Utilities.LFP.clusterStats(ppc_C1, ppc_C2, frqs);
       

%         keyboard
    end
end












