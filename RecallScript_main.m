

%% Recall script
% same vein as AIC script
% need to create a code package to do the same
% Need a setPathsRecall script/function
% Extractcells
% Preprocess data - make rasters, carve up different sections
% Think about analysis - Do recall FR's obey the same rank ordering as encoding?



%% Carving up recall sections

% With current task structure each trial has 2 free recall sections
% Trial Begin
% Show all images
% Visual search distractor
% Free recall part 1
%     Q1: What did you see? Q2: Describe them
% Cued recall (My saying which to imagine)
% Free recall part 2
%     Q3: Can you describe them again

%% Set paths to stuff
setDiskPaths

paths.basePath = diskPath;
% paths.basePath = 'Z:\LabUsers\vwadia\SUAnalysis';

% paths.patientPath = 'P71CS';
% paths.sessPath = 'Recall_Session_1_20201121';
% paths.sessPath = 'Recall_Session_2_20201124';

% paths.patientPath = 'P76CS';
% paths.sessPath = 'ReScreenRecall_Session_1_20210917';
% paths.sessPath = 'Recall_Session_2_20210925';
% paths.sessPath = 'ReScreenRecall_Session_3_20210927';

% paths.patientPath = 'P79CS';
% paths.sessPath = 'ReScreenRecall_Session_1_20220330';
% paths.sessPath = 'ReScreenRecall_Session_2_20220403';
% paths.sessPath = 'ReScreenRecall_Session_3_20220405';

% paths.patientPath = 'P80CS';
% paths.sessPath = 'ReScreenRecall_Session_1_20220728';
% paths.sessPath = 'ReScreenRecall_Session_2_20220731';

% paths.patientPath = 'P84CS';
% paths.sessPath = 'ReScreenRecall_Session_1_20230406';
% paths.sessPath = 'ReScreenRecall_Session_2_20230408';

paths.patientPath = 'P85CS';
paths.sessPath = 'ReScreenRecall_Session_1_20230424';


paths.rawPath = 'raw'; 

paths.sortPath = 'sort';

paths.taskPath = 'Recall_Task';

pathFrames = [paths.basePath filesep paths.taskPath filesep paths.patientPath filesep paths.sessPath filesep 'stimuliUsedRecall'];

basePath = [paths.basePath filesep paths.taskPath filesep paths.patientPath filesep paths.sessPath];

% if strcmp(host(1:end-1), 'DWA644201')
%     paths.taskCodePath = 'D:\Users\wadiav\Dropbox\Caltech\Thesis\Human_work\Cedars\RecallTaskVarun';
% elseif strcmp(host(1:end-1), 'DESKTOP-LJHLIED')
%     paths.taskCodePath = 'E:\Dropbox\Caltech\Thesis\Human_work\Cedars\RecallTaskVarun';
% elseif strcmp(host(1:end-1), 'Varuns-iMac-2.local')
%     paths.taskCodePath = '/Users/varunwadia/Dropbox/Caltech/Thesis/Human_work/Cedars/RecallTaskVarun';
% elseif strcmp(host(1:end-1), 'MacBook-Pro-3.local')
%     paths.taskCodePath = '/Volumes/Macintosh HD/Users/varunwadia/Dropbox/Caltech/Thesis/Human_work/Cedars/RecallTaskVarun';    
% end
paths.taskCodePath = [boxPath filesep 'RecallTaskVarun'];

paths.dataPath = 'processedData';

% addpath(genpath([paths.basePath filesep 'helpers']));
addpath(genpath([paths.basePath filesep 'Code' filesep 'osortTextUI']));
addpath(genpath([paths.basePath filesep 'ObjectSpace']));
addpath(paths.taskCodePath);

% [paths.taskCodePath filesep setTTLCodes];
setTTLCodes


% get rolling already
% datStrct = 'RecallData.mat'; AllCells = 1;
% datStrct = 'RecallData_NoFreeRec.mat'; AllCells = 1;
% datStrct = 'RecallData_FFAonly.mat'; IT_Cells = 1;
% datStrct = 'RecallData_IT_MTL.mat'; IT_MTL_Cells = 1;
% 
% if exist([paths.basePath filesep paths.taskPath filesep paths.patientPath filesep paths.sessPath filesep datStrct], 'file')
%     load([paths.basePath filesep paths.taskPath filesep paths.patientPath filesep paths.sessPath filesep datStrct])
%     strctCELL = struct2cell(strctCells');
%     strctCELL = strctCELL';
% end

%% plot some shit - per cell all stim

% Recall.recall_plotPerCellAllStim
%  
% % plot some more shit - per stim all cells
% 
% Recall.recall_plotPerStimAllCells
% 
% % plot even more shit - per cell per stim
% 
% Recall.recall_plotPerCellPerStim
% 
% % Plotting fulltrialFR - all cells
% 
% Recall.recall_plotFullTrialFRAllCells

%%

FFAChansOnly = 0; 

paths.cellPath = [paths.basePath filesep paths.taskPath filesep paths.patientPath filesep paths.sessPath filesep paths.sortPath filesep 'final'];
events = getRawTTLs([paths.basePath filesep paths.taskPath filesep paths.patientPath filesep paths.sessPath filesep paths.rawPath filesep 'Events.nev'], 1); % TTLs
RecallData.eventsMS = events;
RecallData.eventsMS(:, 1) = events(:, 1)*1e-3;

if strcmp(paths.sessPath,'ReScreenRecall_Session_1_20210917')
    exp_start = find(RecallData.eventsMS(:, 2) == EXPERIMENT_ON);
    exp_end = find(RecallData.eventsMS(:, 2) == EXPERIMENT_OFF);
    RecallData.eventsMS = RecallData.eventsMS(exp_start(2):exp_end(3), :);
elseif strcmp(paths.sessPath,'ReScreenRecall_Session_3_20210927')
    exp_start = find(RecallData.eventsMS(:, 2) == EXPERIMENT_ON);
    exp_end = find(RecallData.eventsMS(:, 2) == EXPERIMENT_OFF);
    RecallData.eventsMS = RecallData.eventsMS(exp_start(1):exp_end(1), :);
elseif strcmp(paths.sessPath,'ReScreenRecall_Session_1_20220330')    
    exp_start = find(RecallData.eventsMS(:, 2) == EXPERIMENT_ON);
    exp_end = find(RecallData.eventsMS(:, 2) == EXPERIMENT_OFF);
    RecallData.eventsMS = RecallData.eventsMS(exp_start(2):exp_end(1), :);  
elseif strcmp(paths.sessPath,'ReScreenRecall_Session_2_20220403')    
    exp_start = find(RecallData.eventsMS(:, 2) == EXPERIMENT_ON);
    exp_end = find(RecallData.eventsMS(:, 2) == EXPERIMENT_OFF);
    RecallData.eventsMS = RecallData.eventsMS(exp_start(1):exp_end(1), :);    
elseif strcmp(paths.sessPath,'ReScreenRecall_Session_3_20220405')    
    exp_start = find(RecallData.eventsMS(:, 2) == EXPERIMENT_ON);
    exp_end = find(RecallData.eventsMS(:, 2) == EXPERIMENT_OFF);
    RecallData.eventsMS = RecallData.eventsMS(exp_start(2):exp_end(2), :);    
elseif strcmp(paths.sessPath,'ReScreenRecall_Session_1_20220728')
    exp_start = find(RecallData.eventsMS(:, 2) == EXPERIMENT_ON);
    exp_end = find(RecallData.eventsMS(:, 2) == EXPERIMENT_OFF);
    RecallData.eventsMS = RecallData.eventsMS(exp_start(1):exp_end(1), :); 
elseif strcmp(paths.sessPath,'ReScreenRecall_Session_2_20220731')
    exp_start = find(RecallData.eventsMS(:, 2) == EXPERIMENT_ON);
    exp_end = find(RecallData.eventsMS(:, 2) == EXPERIMENT_OFF);
    RecallData.eventsMS = RecallData.eventsMS(exp_start(1):exp_end(1), :); 
elseif strcmp(paths.sessPath,'ReScreenRecall_Session_1_20230406') ||...
        strcmp(paths.sessPath,'ReScreenRecall_Session_2_20230408')
    exp_start = find(RecallData.eventsMS(:, 2) == EXPERIMENT_ON);
    exp_end = find(RecallData.eventsMS(:, 2) == EXPERIMENT_OFF);
    RecallData.eventsMS = RecallData.eventsMS(exp_start(1):exp_end(1), :);
elseif strcmp(paths.sessPath,'ReScreenRecall_Session_1_20230424')
    exp_start = find(RecallData.eventsMS(:, 2) == EXPERIMENT_ON);
    exp_end = find(RecallData.eventsMS(:, 2) == EXPERIMENT_OFF);
    RecallData.eventsMS = RecallData.eventsMS(exp_start(3):exp_end(1), :);
else
    exp_start = find(RecallData.eventsMS(:, 2) == EXPERIMENT_ON);
    exp_end = find(RecallData.eventsMS(:, 2) == EXPERIMENT_OFF);
    % do this better - i.e. per session
    if length(exp_start) > 1
        if length(exp_end) > 1
            RecallData.eventsMS = RecallData.eventsMS(exp_start(2):exp_end(2), :);
        else
            RecallData.eventsMS = RecallData.eventsMS(exp_start(2):exp_end, :);
        end
    else
        RecallData.eventsMS = RecallData.eventsMS(exp_start:exp_end, :);
    end
end

if ~exist('strctCells')
    [strctCells, dupCells] = Utilities.extractCells(paths.cellPath, [paths.basePath filesep paths.taskPath filesep paths.patientPath filesep paths.sessPath]);
end

%% carve up different brain regions  -  for state space analyses

% Could also just input the appropriate channels to 'extractCells'

strctCELL = struct2cell(strctCells');
strctCELL = strctCELL';

% IT by itself
% IT+MTL
% IT+MTL+MFC
% IT+MFC

if FFAChansOnly
    IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
    strctCells = strctCells(IT_Cells);
end

% IT_Cells = cellfun(@(x) strcmp(x, 'RFFA'), strctCELL(:, 4));
% strctCells = strctCells(IT_Cells);

% IT_MTL_Cells = cellfun(@(x) ismember(x, {'LA', 'LH', 'RA', 'RH', 'RFFA'}), strctCELL(:, 4), 'UniformOutput', false);
% strctCells = strctCells(cell2mat(IT_MTL_Cells));

% IT_MTL_MFC_Cells = cellfun(@(x) strcmp(x, {'LA', 'LH', 'RA', 'RH', 'LSMA', 'RSMA', 'RFFA'}), strctCELL(:, 4));
% strctCells = strctCells(cell2mat(IT_MTL_MFC_Cells));

% IT_MFC_Cells = cellfun(@(x) strcmp(x, ['LSMA', 'RSMA', 'RFFA']), strctCELL(:, 4));
% strctCells = strctCells(cell2mat(IT_MFC_Cells));

%% Best way to arrange data

% make each trial a standalone
% EACH image needs its own TTL. Write the image name to a log file (to have
% ground truth)
train_start = find(RecallData.eventsMS(:, 2) == TRAINING_RECALL_BEGIN); 
train_end = find(RecallData.eventsMS(:, 2) == TRAINING_RECALL_END);
RecallData.eventsMS = RecallData.eventsMS(train_end:end, :);

RecallData = Recall.preProcessRecall(strctCells, paths, RecallData, atCedars);

RecallData.offsetTones = [1000 6500]; % no stimOFFtimes used so raster is -offset1 to +offset2 around the event
RecallData.offsetFR = [5000 3000]; % no stimOFFtimes used so raster is -offset1 to +offset2 around the event
RecallData.offsetEnc = [1500 500];

% number of FR repetitions
% RecallData.FR_numRepetitions = length(RecallData.recEventMarkers{1, 3}) + length(RecallData.recEventMarkers{1, 5}); 
RecallData.CR_numRepetitions = length(RecallData.CROrder)/length(unique(RecallData.CROrder));
RecallData.Enc_numRepetitions = length(RecallData.EncodingOrder)/length(unique(RecallData.EncodingOrder));
%%
bin_size = 200; % ms - for smoothing
use_both_offsets = 1;
order = 1;
% carve up rasters for encoding, CR and FR - 496s for P71 Session 2
% carve up rasters for encoding, CR and FR - 800s for P76 Session 1 - went by in 27 seconds on Feb4th...is this wrong? (nope only FFA cells)
% carve up rasters for encoding, CR and FR - s for P76 Session 2
% carve up rasters for encoding, CR and FR - s for P76 Session 3
tic
for cellIndex = 1:length(strctCells)
    
%     for stim = 1:length(RecallData.stimuli) % there are 2 recall sections per trial        
%         [RecallData.FREventTimeCourse{stim, 1}{cellIndex, 1}, RecallData.FREventTimeCourse{stim, 1}{cellIndex, 2},  RecallData.FREventTimeCourse{stim, 1}{cellIndex, 3}]...
%             = Utilities.makeRastersFromStrctCells(strctCells(cellIndex), RecallData.offsetFR, use_both_offsets, order, bin_size, cell2mat(RecallData.recEventMarkers{stim, 6}));       
%     end
%     
%     for trl = 1:length(RecallData.trialONTimes) % read in transcript for each trial
%         [RecallData.FRTimeCourse{trl, 1}{cellIndex, 1}, RecallData.FRTimeCourse{trl, 1}{cellIndex, 2},  RecallData.FRTimeCourse{trl, 1}{cellIndex, 3}]...
%             = Utilities.makeRastersFromStrctCells(strctCells(cellIndex), RecallData.offsetFR, use_both_offsets, order, bin_size, RecallData.audioStartTimes(trl), RecallData.audioEndTimes(trl));        
%     end
    
    [RecallData.CRTimeCourse{cellIndex, 1}, RecallData.CRTimeCourse{cellIndex, 2},  RecallData.CRTimeCourse{cellIndex, 3}]...
        = Utilities.makeRastersFromStrctCells(strctCells(cellIndex), RecallData.offsetTones, use_both_offsets, RecallData.CROrder, bin_size, RecallData.toneONTimes);
    
    [RecallData.EncodingTimeCourse{cellIndex, 1}, RecallData.EncodingTimeCourse{cellIndex, 2},  RecallData.EncodingTimeCourse{cellIndex, 3}]...
        = Utilities.makeRastersFromStrctCells(strctCells(cellIndex), RecallData.offsetEnc, 1, 1, bin_size, RecallData.imageONTimes, RecallData.imageOFFTimes);
    
    if ~strcmp(paths.sessPath, 'ReScreenRecall_Session_2_20230408') && ~strcmp(paths.sessPath, 'ReScreenRecall_Session_1_20230424') % these sessions were done without distraction period
        [RecallData.DistractionTimeCourse{cellIndex, 1}, RecallData.DistractionTimeCourse{cellIndex, 2},  RecallData.DistractionTimeCourse{cellIndex, 3}]...
            = Utilities.makeRastersFromStrctCells(strctCells(cellIndex), RecallData.offsetEnc, 1, 1, bin_size, RecallData.distONTimes, RecallData.distOFFTimes);
    end
    
    % for reactivation calculation
     [RecallData.PreTrialBaselineTimeCourse{cellIndex, 1}, RecallData.PreTrialBaselineTimeCourse{cellIndex, 2},  RecallData.PreTrialBaselineTimeCourse{cellIndex, 3}]...
        = Utilities.makeRastersFromStrctCells(strctCells(cellIndex), RecallData.offsetFR, 1, 1, bin_size, RecallData.trialONTimes);
    
       % for reactivation calculation - need to assign the ton-period-begin time in the preprocess script
     [RecallData.PreCRBaselineTimeCourse{cellIndex, 1}, RecallData.PreCRBaselineTimeCourse{cellIndex, 2},  RecallData.PreCRBaselineTimeCourse{cellIndex, 3}]...
        = Utilities.makeRastersFromStrctCells(strctCells(cellIndex), [5000 5000], 1, 1, bin_size, RecallData.tonePeriodBeginTimes);
end
toc
%%

if ~strcmp(paths.sessPath, 'Recall_Session_2_20210925')
    Recall.recall_plotPerCellRasterandBarPerStim
elseif strcmp(paths.sessPath, 'Recall_Session_2_20210925') && FFAChansOnly
    Recall.recall_plotPerCellRasterandBarPerStim
end
keyboard % safety in case it doesn't enter script for whatever reason

% go to recall_plotPerCellRasterandBarPerStim for making responses per
% image
% save([paths.basePath filesep paths.taskPath filesep paths.patientPath filesep paths.sessPath filesep 'RecallData_NoFReeRec.mat'], 'RecallData', 'strctCells', '-v7.3')

%% collecting population responses - free recall

 RecallData.perCellAllStimFR = cell(length(strctCells), 3);
 RecallData.order_perCellAllStimFR = repelem(1:length(RecallData.stimuli), RecallData.FR_numRepetitions)';
% per cell all stim
for rep = 1:length(strctCells)
    for stim = 1:length(RecallData.stimuli)
        RecallData.perCellAllStimFR{rep, 1} = [RecallData.perCellAllStimFR{rep, 1}; RecallData.FREventTimeCourse{stim, 1}{rep, 1}(:, 1:sum(RecallData.offsetFR))];
        RecallData.perCellAllStimFR{rep, 2} = [RecallData.perCellAllStimFR{rep, 2}; RecallData.FREventTimeCourse{stim, 1}{rep, 2}(:, 1:sum(RecallData.offsetFR))];
        if stim == 1 
            RecallData.perCellAllStimFR{rep, 3} = [RecallData.perCellAllStimFR{rep, 3}; RecallData.FREventTimeCourse{stim, 1}{rep, 3}(:, 1:sum(RecallData.offsetFR))];
        end
    end
end

RecallData.perStimAllCellsFR = cell(length(RecallData.stimuli), 3);
RecallData.order_perStimAllCellsFR = repelem(1:length(strctCells), RecallData.FR_numRepetitions)';

% per stim all cells
for stim = 1:length(RecallData.stimuli)
    for rep = 1:length(strctCells)
        RecallData.perStimAllCellsFR{stim, 1} = [RecallData.perStimAllCellsFR{stim, 1}; RecallData.FREventTimeCourse{stim, 1}{rep, 1}(:, 1:sum(RecallData.offsetFR))];
        RecallData.perStimAllCellsFR{stim, 2} = [RecallData.perStimAllCellsFR{stim, 2}; RecallData.FREventTimeCourse{stim, 1}{rep, 2}(:, 1:sum(RecallData.offsetFR))];
        if rep == 1
            RecallData.perStimAllCellsFR{stim, 3} = [RecallData.perStimAllCellsFR{stim, 3}; RecallData.FREventTimeCourse{stim, 1}{rep, 3}(:, 1:sum(RecallData.offsetFR))];
            
        end
    end
end
% collect population responses 

% CR
RecallData.perStimAllCellsCR = cell(length(RecallData.stimuli), 3);
RecallData.order_perStimAllCellsCR = repelem(1:length(strctCells), RecallData.CR_numRepetitions)';

% per stim all cells
for stim = 1:length(RecallData.stimuli)
    for rep = 1:length(strctCells)
        stimRaster = RecallData.CRTimeCourse{rep, 1}(find(RecallData.CROrder == stim), 1:sum(RecallData.offsetTones));
        stimPsth = RecallData.CRTimeCourse{rep, 2}(find(RecallData.CROrder == stim), 1:sum(RecallData.offsetTones));
        times = RecallData.CRTimeCourse{rep, 3}; 
        RecallData.perStimAllCellsCR{stim, 1} = [RecallData.perStimAllCellsCR{stim, 1}; stimRaster];
        RecallData.perStimAllCellsCR{stim, 2} = [RecallData.perStimAllCellsCR{stim, 2}; stimPsth];
        
        if rep == 1
            RecallData.perStimAllCellsCR{stim, 3} = [RecallData.perStimAllCellsCR{stim, 3}; times];    
        end
    end
end
%
% Encoding
RecallData.perStimAllCellsEncoding = cell(length(RecallData.stimuli), 3);
RecallData.order_perStimAllCellsEncoding = repelem(1:length(strctCells), RecallData.Enc_numRepetitions)';

% per stim all cells
for stim = 1:length(RecallData.stimuli)
    for rep = 1:length(strctCells)
        stimRaster = RecallData.EncodingTimeCourse{rep, 1}(find(RecallData.EncodingOrder == stim), :);
        stimPsth = RecallData.EncodingTimeCourse{rep, 2}(find(RecallData.EncodingOrder == stim), :);
        times = RecallData.EncodingTimeCourse{rep, 3}; 
        if rep == 11 && strcmp(paths.sessPath, 'ReScreenRecall_Session_3_20210927')
            RecallData.perStimAllCellsEncoding{stim, 1} = [RecallData.perStimAllCellsEncoding{stim, 1}; stimRaster(:, 1:end-1)];
            RecallData.perStimAllCellsEncoding{stim, 2} = [RecallData.perStimAllCellsEncoding{stim, 2}; stimPsth(:, 1:end-1)];
        else
            RecallData.perStimAllCellsEncoding{stim, 1} = [RecallData.perStimAllCellsEncoding{stim, 1}; stimRaster];
            RecallData.perStimAllCellsEncoding{stim, 2} = [RecallData.perStimAllCellsEncoding{stim, 2}; stimPsth];
        end
        
        if rep == 1
            RecallData.perStimAllCellsEncoding{stim, 3} = [RecallData.perStimAllCellsEncoding{stim, 3}; times];    
        end
    end
end



%% Full trial FR
strctCELL = struct2cell(strctCells');
strctCELL = strctCELL';
RecallData.perTrialAllCellsFullFR = cell(length(RecallData.trialONTimes), 3);
RecallData.order_perTrialAllCellsFullFR = strctCELL(:, 3);

% per stim all cells
for trl = 1:length(RecallData.trialONTimes)
    for cellIndex = l(strctCells)
%         if trl == 3
%             RecallData.perTrialAllCellsFullFR{trl, 1} = [RecallData.perTrialAllCellsFullFR{trl, 1}; RecallData.FRTimeCourse{trl, 1}{cellIndex, 1}(1:98073)];
%             RecallData.perTrialAllCellsFullFR{trl, 2} = [RecallData.perTrialAllCellsFullFR{trl, 2}; RecallData.FRTimeCourse{trl, 1}{cellIndex, 2}(1:98073)];
%         else
            RecallData.perTrialAllCellsFullFR{trl, 1} = [RecallData.perTrialAllCellsFullFR{trl, 1}; RecallData.FRTimeCourse{trl, 1}{cellIndex, 1}];
            RecallData.perTrialAllCellsFullFR{trl, 2} = [RecallData.perTrialAllCellsFullFR{trl, 2}; RecallData.FRTimeCourse{trl, 1}{cellIndex, 2}];
%         end
        if cellIndex == 1
            RecallData.perTrialAllCellsFullFR{trl, 3} = [RecallData.perTrialAllCellsFullFR{trl, 3}; RecallData.FRTimeCourse{trl, 1}{cellIndex, 3}];
            
        end
    end
end

keyboard

save([paths.basePath filesep paths.taskPath filesep paths.patientPath filesep paths.sessPath filesep 'RecallData.mat'], 'RecallData', 'strctCells', '-v7.3')


%% carving up audio per trial
%% commented out after single use
% to collect audio simply use the frist trial recall onset and the
% experiment end as markers
%{
TRIALON = find(RecallData.eventsMS(:, 2) == TRIAL_ON);
numTrials = length(TRIALON);

recallON = find(RecallData.eventsMS(:, 2) == DISTRACTION_PERIOD_END);
recallOFF = [TRIALON(2:end); find(RecallData.eventsMS(:, 2) == EXPERIMENT_OFF)];
recallOFFTimes = RecallData.eventsMS(recallOFF, 1); % resetting these here
recallONTimes = RecallData.eventsMS(recallON, 1);

audioPath = [paths.basePath filesep paths.taskPath filesep paths.patientPath filesep paths.sessPath];
Fs = 32000; % Nlx micro sampling rate
audioFileMS = Utilities.makeAudioFile([audioPath filesep paths.rawPath], 190);
% audioFileMS(:, 2) = highpass(audioFileMS(:, 2), 250, Fs);
for trl = 1:numTrials
    
    audioInRange = and(audioFileMS(:, 1) >= recallONTimes(trl), audioFileMS(:, 1) <= recallOFFTimes(trl));
    recallSection = audioFileMS(audioInRange, 2)/Fs;
   
    suff = ['_' num2str(trl) '.wav']; 
    audiowrite([audioPath filesep 'FreeRecallSection' suff], recallSection, Fs);
    % to play
    % soundsc(recallSection, Fs);
    
end
%}