
%%

% Run RecallScript until first breakpoint (making CR and Enc time courses
% Then come here
% April 2022

%% compute response latency of cell for encoding - could just take this from screeningScript

% collect images
imDir = dir(fullfile(pathFrames));
imDir = imDir(~ismember({imDir.name}, {'.', '..', '.DS_Store', 'Thumbs.db'}));

% making it same form as screening script
RecallData.imageIDs = unique(RecallData.EncodingOrder);
% can get this from psth{cellIndex, 3} or the 'times' variable
RecallData.timelimits = [-RecallData.offsetEnc(1) (size(RecallData.EncodingTimeCourse{1, 1}, 2)-RecallData.offsetEnc(1))]*1e-3;

imageOffPoints = find(RecallData.eventsMS(:, 2) == IMAGE_OFF);
imageOnPoints = find(RecallData.eventsMS(:, 2) == IMAGE_ON);

imageOffTimes = RecallData.eventsMS(imageOffPoints, 1);
imageOnTimes = RecallData.eventsMS(imageOnPoints, 1);

if isequal(length(imageOffTimes), length(imageOnTimes))
    stimDur = median(imageOffTimes - imageOnTimes)*1e-3;
    stimOffDur = (imageOnTimes(2:end) - imageOffTimes(1:end-1));
    stimOffDur = median(stimOffDur(stimOffDur < 1e6))*1e-3;
else
    stimDur = (imageOffTimes(1) - imageOnTimes(1))*1e-3;
    stimOffDur = (imageOnTimes(2) - imageOffTimes(1))*1e-3;
end
RecallData.stimDur = stimDur;
RecallData.stimOffDur = stimOffDur;

%% collecting screening spike counts for relevant images
AllCells = 0;
IT_MTL_Cells_Only = 0;

% P76
if strcmp(paths.sessPath, 'ReScreenRecall_Session_1_20210917')
    options.recalled_stim = [12 19 25 123 270 487];
    
elseif strcmp(paths.sessPath, 'Recall_Session_2_20210925') % have to handle this separately
    options.recalled_stim = [54 129 130 186 270 449];
    
elseif strcmp(paths.sessPath, 'ReScreenRecall_Session_3_20210927')
    options.recalled_stim = [18 44 45 81 135 181 230 344] ;
 % P79   
elseif strcmp(paths.sessPath, 'ReScreenRecall_Session_1_20220330')
    options.recalled_stim = [9 157 167 200 201 291 422 498];
    
elseif strcmp(paths.sessPath, 'ReScreenRecall_Session_2_20220403')
    options.recalled_stim = [9 12 117 292 360 368 421 492];
    
elseif strcmp(paths.sessPath, 'ReScreenRecall_Session_3_20220405')
    options.recalled_stim = [77 112 160 232 278 345 387 440];
  % P80  
elseif strcmp(paths.sessPath, 'ReScreenRecall_Session_1_20220728')
    options.recalled_stim = [17 61 76 114 157 161 177 480 ];

elseif strcmp(paths.sessPath, 'ReScreenRecall_Session_2_20220731')
    options.recalled_stim = [55 88 148 251 256 274 285 365];
    
elseif strcmp(paths.sessPath, 'ReScreenRecall_Session_1_20230406')
    options.recalled_stim = [68 121 243 261 281 308 415 434];

elseif strcmp(paths.sessPath, 'ReScreenRecall_Session_2_20230408')
    options.recalled_stim = [52 141 195 223 238 351 399 482];
    
elseif strcmp(paths.sessPath, 'ReScreenRecall_Session_1_20230424')
    options.recalled_stim = [7 11 99 153 201 251 355 486];
end

if strcmp(paths.sessPath, 'Recall_Session_2_20210925')
    load([diskPath filesep paths.taskPath  filesep paths.patientPath filesep 'RecallScreening_Session_2_20210925' filesep 'PsthandResponses']);
    
    if FFAChansOnly
        s_cells = load([diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'RecallScreening_Session_2_20210925' filesep 'strctCells']);
        s_CELL = struct2cell(s_cells.strctCells');
        s_CELL = s_CELL';
        
        IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), s_CELL(:, 4));
        psths = psths(IT_Cells, :);
        responses = responses(IT_Cells, :);
    end
    
    
    pathOut = [basePath filesep 'processedData' filesep 'SpikeCountComparison_MorningScreenandCR'];
    
else
    load([paths.basePath filesep paths.taskPath filesep paths.patientPath filesep paths.sessPath filesep 'PsthandResponses']);
    pathOut = [basePath filesep 'processedData' filesep 'SpikeCountComparison_ReScreenandCR'];
    
end

strctCELL = struct2cell(strctCells');
strctCELL = strctCELL';


if FFAChansOnly && ~strcmp(paths.sessPath, 'Recall_Session_2_20210925')
    psths = psths(IT_Cells, :);
    responses = responses(IT_Cells, :);
elseif IT_MTL_Cells_Only
    IT_MTL_Cells = cellfun(@(x) ismember(x, {'LA', 'LH', 'RA', 'RH', 'RFFA', 'LFFA'}), strctCELL(:, 4), 'UniformOutput', false);
    strctCells = strctCells(cell2mat(IT_MTL_Cells));
    psths = psths(IT_MTL_Cells, :);
    responses = responses(IT_MTL_Cells, :);
end


% old way
% index = cellfun(@isempty, responses(:, 2));

% GETTING RID OF THIS ----------------------------------------------------
% index = cell2mat(cellfun(@isnan, responses(:, 2), 'UniformOutput', false));
% responses(index(:, 1), :) = []; % remove non-responsive neurons
% if ~strcmp(paths.sessPath, 'Recall_Session_2_20210925') % for this session can only do IT neurons as the rest are not matchable
%     strctCells(index(:, 1)) = [];
%     RecallData.CRTimeCourse(index(:, 1), :) = [];
%     RecallData.EncodingTimeCourse(index(:, 1), :) = [];
%     RecallData.DistractionTimeCourse(index(:, 1), :) = [];
% end
% -------------------------------------------------------------------------
RecallData.ScreenResponses = responses;

% grab the screening responses for only the imagined stimuli
for cellIndex = 1:length(strctCells)
    if ~isnan(responses{cellIndex, 1})
        RecallData.ScreenResponses{cellIndex, 4} = responses{cellIndex, 1}(options.recalled_stim, 1);
    else
        RecallData.ScreenResponses{cellIndex, 4} = nan;
    end
end

%% calculating the 'encoding' responses (not to be confused with the 'screening' spike counts)
if strcmp(paths.sessPath, 'Recall_Session_2_20210925')
    method = 0; % for computing from encoding use method 0
else
    method = 3; % else from screening poisson
end
for cellIndex = 1:length(strctCells)
    
    % for encoding we can grab this from screening data
    if strcmp(paths.sessPath, 'Recall_Session_2_20210925')
        switch method
            case 0
                n_stdDevs = 2.5;
                [respLat, max_group] = Utilities.computeResponseLatency(RecallData.EncodingTimeCourse(cellIndex, :), RecallData.EncodingOrder, RecallData.timelimits,...
                    RecallData.stimOffDur, RecallData.stimDur, method, n_stdDevs);
                RecallData.EncResponses{cellIndex, 2} = respLat - (-RecallData.timelimits(1)*1e3);
                
            case 1
                [respLat, max_group] = Utilities.computeResponseLatency(RecallData.EncodingTimeCourse(cellIndex, :), RecallData.EncodingOrder, RecallData.timelimits,...
                    RecallData.stimOffDur, RecallData.stimDur);
                RecallData.EncResponses{cellIndex, 2} = respLat - (-RecallData.timelimits(1)*1e3);
                
            case 2
                [respLat, max_group] = Utilities.computeResponseLatency(RecallData.EncodingTimeCourse(cellIndex, :), RecallData.EncodingOrder, RecallData.timelimits,...
                    RecallData.stimOffDur, RecallData.stimDur);
                RecallData.EncResponses{cellIndex, 2} = respLat - (-RecallData.timelimits(1)*1e3);
                
            case 3
                [respLat, ~] = Utilities.computeRespLatPoisson(RecallData.EncodingTimeCourse(cellIndex, :), RecallData.EncodingOrder, RecallData.EncodingOrder, RecallData.timelimits,...
                    RecallData.stimDur, false, 'trial', 2.5);
                RecallData.EncResponses{cellIndex, 2} = respLat; % poisson spits it out already adjusted
                
        end
    else
        respLat = RecallData.ScreenResponses{cellIndex, 2};
        RecallData.EncResponses{cellIndex, 2} = respLat;
    end
    RecallData.EncResponses{cellIndex, 1} = [];

    endRas = size(RecallData.EncodingTimeCourse{cellIndex, 1}, 2);
    if respLat ~= 0 && ~isnan(respLat)
        windowLength = floor(RecallData.stimDur*1e3);
        windowBegin = floor(respLat) + (-RecallData.timelimits(1)*1e3); % this was unadjusted the whole time...bruh WTF. Luckily probably won't change anything vwadia March2023
        windowEnd = windowBegin+windowLength;
        if windowEnd > endRas
            windowEnd = endRas;
        end
        
        for i = l(RecallData.imageIDs)
            EncStimRaster = RecallData.EncodingTimeCourse{cellIndex, 1}(find(RecallData.EncodingOrder == RecallData.imageIDs(i)), windowBegin:windowEnd);
            RecallData.EncResponses{cellIndex, 1}(i, 1) = mean(mean(EncStimRaster))*1e3;
        end
        
        RecallData.EncResponses{cellIndex, 3} = strctCells(cellIndex).Name;
        if exist('max_group')
             RecallData.EncResponses{cellIndex, 4} = max_group;
        end
        respLat = 0; % resetting - don't need to do this really
    end
end
%% computing the CR spike counts for the images
CRLength = 5000; %ms
for cellIndex = 1:length(strctCells)
    spikeCount = [];
    stimRaster = {};
    for stim = 1:length(RecallData.imageIDs)
        full_CR_psth = RecallData.CRTimeCourse(cellIndex, 1:3);
        CR_psth = {full_CR_psth{1, 1}(find(RecallData.CROrder == stim), :), full_CR_psth{1, 2}(find(RecallData.CROrder == stim), :), full_CR_psth{1, 3}};
        CRStimRaster = CR_psth{1, 1};
        stimRaster{stim, 1} = CRStimRaster;
        spikeCount(end+1, :) = mean(mean(CRStimRaster(:, RecallData.offsetTones:RecallData.offsetTones+CRLength)))*1e3;
    end
    RecallData.CRResponses{cellIndex, 1} = spikeCount;
%     test{cellIndex, 1} = spikeCount;
end
    

strctResp = struct; 
for cellIndex = 1:length(strctCells)
    strctResp(cellIndex).Name = strctCells(cellIndex).Name;
    strctResp(cellIndex).CRResp =  RecallData.CRResponses{cellIndex, 1};
    strctResp(cellIndex).EncResp =  RecallData.EncResponses{cellIndex, 1};
    strctResp(cellIndex).ScrnResp =  RecallData.ScreenResponses{cellIndex, 4};
    
end

% keyboard
if FFAChansOnly
    save([basePath filesep 'ITResponses'],'strctResp');
elseif IT_MTL_Cells_Only 
    save([basePath filesep 'IT_MTLResponses'],'strctResp');
else
    save([basePath filesep 'AllResponses'],'strctResp');   
end

% was saving full recall data (unless only IT neurons)
if ~FFAChansOnly
    save([paths.basePath filesep paths.taskPath filesep paths.patientPath filesep paths.sessPath filesep 'RecallData_NoFReeRec.mat'], 'RecallData', 'strctCells', '-v7.3')
    disp("Saved!")
end

%% old way - Don't really do this anyumore MArch2023
%{
keyboard
for cellIndex = 1:length(strctCells)
   RecallData.CRResponses{cellIndex, 1} = strctResp(cellIndex).CRResp;
   RecallData.EncResponses{cellIndex, 1} = strctResp(cellIndex).EncResp;
   RecallData.ScreenResponses{cellIndex, 4} = strctResp(cellIndex).ScrnResp;
    
end
% CRResponses = RecallData.CRResponses;
% ScrnResponses = RecallData.ScreenResponses(:, 4);
% EncResponses = RecallData.EncResponses;
% save([basePath filesep 'AllResponses'], 'CRResponses', 'ScrnResponses', 'EncResponses');
%% Plot encoding vs recall on the same plot


for cellIndex = l(strctCells)
    spikeCount = {};
    for stim = 1:length(RecallData.stimuli)
        
        % collect psths
        full_enc_psth = RecallData.EncodingTimeCourse(cellIndex, 1:3);
        full_CR_psth = RecallData.CRTimeCourse(cellIndex, 1:3);
        
        enc_psth = {full_enc_psth{1, 1}(find(RecallData.EncodingOrder == stim), :), full_enc_psth{1, 2}(find(RecallData.EncodingOrder == stim), :), full_enc_psth{1, 3}};
        CR_psth = {full_CR_psth{1, 1}(find(RecallData.CROrder == stim), :), full_CR_psth{1, 2}(find(RecallData.CROrder == stim), :), full_CR_psth{1, 3}};
        
        
        EncStimRaster = enc_psth{1, 1};
        CRStimRaster = CR_psth{1, 1};
        stimRaster = {EncStimRaster, CRStimRaster};
        spikeCount(end+1, :) = {RecallData.EncResponses{cellIndex, 1}(stim), mean(mean(CRStimRaster(:, RecallData.offsetTones:end)))*1e3};


    end
        
    
    
    for stim = 1:length(RecallData.imageIDs)
        % compute stimraster
%         EncStimRaster = RecallData.EncodingTimeCourse{cellIndex, 1}(find(RecallData.EncodingOrder == RecallData.imageIDs(stim)), :);
        % collect psths
        % this repetition is stupid find a better way
        % -----------------------------------------------------------------
        full_enc_psth = RecallData.EncodingTimeCourse(cellIndex, 1:3);
        full_CR_psth = RecallData.CRTimeCourse(cellIndex, 1:3);
%         full_FR_psth = RecallData.perCellAllStimFR(cellIndex, 1:3);
        
        enc_psth = {full_enc_psth{1, 1}(find(RecallData.EncodingOrder == stim), :), full_enc_psth{1, 2}(find(RecallData.EncodingOrder == stim), :), full_enc_psth{1, 3}};
        CR_psth = {full_CR_psth{1, 1}(find(RecallData.CROrder == stim), :), full_CR_psth{1, 2}(find(RecallData.CROrder == stim), :), full_CR_psth{1, 3}};
%         FR_psth = {full_FR_psth{1, 1}(find(RecallData.order_perCellAllStimFR == stim), :), full_FR_psth{1, 2}(find(RecallData.order_perCellAllStimFR == stim), :), full_FR_psth{1, 3}};
        
        EncStimRaster = enc_psth{1, 1};
        CRStimRaster = CR_psth{1, 1};
        stimRaster = {EncStimRaster, CRStimRaster};
        %------------------------------------------------------------------

        % relevant variables
        imPath = [imDir(RecallData.imageIDs(stim)).folder filesep imDir(RecallData.imageIDs(stim)).name];
        pathOut = [basePath filesep 'processedData' filesep 'SpikeCountComparison'];
        if ~exist(pathOut)
            mkdir(pathOut);
        end
        filename = [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber)...
            '_' num2str(strctCells(cellIndex).Name) '_Stim_' num2str(stim)];
        plotOptions.EncTimelimits = RecallData.timelimits;
        plotOptions.CRTimelimits = [-RecallData.offsetTones(1) RecallData.offsetTones(2)]*1e-3;
        plotOptions.globalyl = max(cell2mat(spikeCount(:)));
%         plotOptions.task = 'Recall';
        
        % produce figure
        handlesToFig = Utilities.Plotting.PlotRasterAndBar_TwoStim(stimRaster, spikeCount(stim, :), imPath, plotOptions);
        sgtitle({[strctCells(cellIndex).brainArea ' ' num2str(strctCells(cellIndex).Name)],...
            ['Spike Count Comparison Stim ' num2str(stim)]})
%         keyboard
        % save
        print(handlesToFig, filename, '-dpng', '-r0');
        close all
    end
     
end

%% plot bar plots for encoding and recall separately

% condition = 'Encoding';
condition = 'CR';
% condition = 'Screening';

for cellIndex = l(strctCells)
    spikeCount = [];
    stimRaster = {};
    for stim = 1:length(RecallData.stimuli)        
        % collect psths
        switch condition
            case 'Encoding'
                full_enc_psth = RecallData.EncodingTimeCourse(cellIndex, 1:3);
                enc_psth = {full_enc_psth{1, 1}(find(RecallData.EncodingOrder == stim), :), full_enc_psth{1, 2}(find(RecallData.EncodingOrder == stim), :), full_enc_psth{1, 3}};
                EncStimRaster = enc_psth{1, 1};
                stimRaster{stim, 1} = EncStimRaster;
                spikeCount(end+1, :) = RecallData.EncResponses{cellIndex, 1}(stim);
                plotOptions.timelimits = RecallData.timelimits;
                plotOptions.color = [0 0.4470 0.7410];
                
            case 'CR'
                full_CR_psth = RecallData.CRTimeCourse(cellIndex, 1:3);               
                CR_psth = {full_CR_psth{1, 1}(find(RecallData.CROrder == stim), :), full_CR_psth{1, 2}(find(RecallData.CROrder == stim), :), full_CR_psth{1, 3}};
                CRStimRaster = CR_psth{1, 1};
                stimRaster{stim, 1} = CRStimRaster;
                spikeCount(end+1, :) = mean(mean(CRStimRaster(:, RecallData.offsetTones:end)))*1e3;
                plotOptions.timelimits = [-RecallData.offsetTones(1) RecallData.offsetTones(2)]*1e-3;
                plotOptions.color = [0.6350 0.0780 0.1840];
            
            case 'Screening'
                full_Screen_psth = psths(cellIndex, 1:3);               
                Screen_psth = {full_Screen_psth{1, 1}(find(order == options.recalled_stim(stim)), :), full_Screen_psth{1, 2}(find(order == options.recalled_stim(stim)), :), full_Screen_psth{1, 3}};
                ScreenStimRaster = Screen_psth{1, 1};
                stimRaster{stim, 1} = ScreenStimRaster;
                spikeCount(end+1, :) = RecallData.ScreenResponses{cellIndex, 4}(stim);
                plotOptions.timelimits = [Screen_psth{1, 3}(1) Screen_psth{1, 3}(end)];
                plotOptions.color = [0.3010 0.7450 0.9330];
                
        end
    end
        
    
    
    for stim = 1:length(RecallData.imageIDs)
        % compute stimraster
        % this repetition is only because of GlobalYl...
        %------------------------------------------------------------------
%         switch condition
%             case 'Encoding'
%                 full_enc_psth = RecallData.EncodingTimeCourse(cellIndex, 1:3);
%                 enc_psth = {full_enc_psth{1, 1}(find(RecallData.EncodingOrder == stim), :), full_enc_psth{1, 2}(find(RecallData.EncodingOrder == stim), :), full_enc_psth{1, 3}};
%                 EncStimRaster = enc_psth{1, 1};
%                 stimRaster = EncStimRaster;
%                 spikeCount(end+1, :) = RecallData.EncResponses{cellIndex, 1}(stim);
%                 plotOptions.timelimits = RecallData.timelimits;
% 
%                 
%             case 'CR'
%                 full_CR_psth = RecallData.CRTimeCourse(cellIndex, 1:3);               
%                 CR_psth = {full_CR_psth{1, 1}(find(RecallData.CROrder == stim), :), full_CR_psth{1, 2}(find(RecallData.CROrder == stim), :), full_CR_psth{1, 3}};
%                 CRStimRaster = CR_psth{1, 1};
%                 stimRaster = CRStimRaster;
%                 spikeCount(end+1, :) = mean(mean(CRStimRaster(:, RecallData.offsetTones:end)))*1e3;
%                 plotOptions.timelimits = [-RecallData.offsetTones(1) RecallData.offsetTones(2)]*1e-3;
% 
%         end
        %------------------------------------------------------------------

        % relevant variables
        imPath = [imDir(RecallData.imageIDs(stim)).folder filesep imDir(RecallData.imageIDs(stim)).name];
        pathOut = [basePath filesep 'processedData' filesep 'SingleStimBarplots'];
        if ~exist(pathOut)
            mkdir(pathOut);
        end
        
        filename = [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber)...
            '_' num2str(strctCells(cellIndex).Name) '_Stim_' num2str(stim) '_' condition];
        plotOptions.globalyl = max(spikeCount(:));
        
        % produce figure
        handlesToFig = Utilities.Plotting.PlotRasterAndBar_SingleStim(stimRaster{stim, 1}, spikeCount(stim, :), imPath, plotOptions);
        sgtitle({[strctCells(cellIndex).brainArea ' ' num2str(strctCells(cellIndex).Name)],...
            ['Spike Count for Stim ' num2str(stim)],...
            [condition]})
        % save
        print(handlesToFig, filename, '-dpng', '-r0');
        close all
    end
     
end
%}
