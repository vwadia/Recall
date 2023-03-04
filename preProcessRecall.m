function [RecallData] = preProcessRecall(strctCells, paths, RecallData, atCedars)
% This function takes in the recallData with some basic things (events etc.),
% which session it is and spits out a full recallData struct with transcripts
% for each trial, block order for images, block order for cued recall etc.
% vwadia June 2021

setTTLCodes;

%% Generic stuff that is session agnostic
% do I need these? - maybe to track cell time course across a whole
% trial?
trialONTimes = RecallData.eventsMS(find(RecallData.eventsMS(:, 2) == TRIAL_ON), 1);
trialOFFTimes = RecallData.eventsMS(find(RecallData.eventsMS(:, 2) == TRIAL_OFF), 1);

imageONTimes = RecallData.eventsMS(find(RecallData.eventsMS(:, 2) == IMAGE_ON), 1);
imageOFFTimes = RecallData.eventsMS(find(RecallData.eventsMS(:, 2) == IMAGE_OFF), 1);

% free recall times - 2 per trial
recallONTimes = RecallData.eventsMS(find(RecallData.eventsMS(:, 2) == RECALL_BEGIN), 1);
recallOFFTimes = sort([trialONTimes(2:end); RecallData.eventsMS(find(RecallData.eventsMS(:, 2) == TONE_PERIOD_BEGIN), 1);...
    RecallData.eventsMS(find(RecallData.eventsMS(:, 2) == EXPERIMENT_OFF), 1)]); % jugaad

% 6 per trial
toneONTimes = sort([RecallData.eventsMS(find(RecallData.eventsMS(:, 2) == TONE_1), 1); RecallData.eventsMS(find(RecallData.eventsMS(:, 2) == TONE_2), 1)]);
toneOFFTimes = toneONTimes + 5000; %ms

distractionONTimes = RecallData.eventsMS(find(RecallData.eventsMS(:, 2) == DISTRACTION_PERIOD_BEGIN), 1);
distractionOFFTimes = RecallData.eventsMS(find(RecallData.eventsMS(:, 2) == DISTRACTION_PERIOD_END), 1);


%% read in transcript and manually choose recall sections
RecallData.FRFULLTS = cell(length(trialONTimes), 1);
RecallData.recEventMarkers = cell(length(trialONTimes), 1);
% How the audio sections were carved up
audioStartTimes = RecallData.eventsMS(find(RecallData.eventsMS(:, 2) == DISTRACTION_PERIOD_END), 1);
audioEndTimes = [];
% for trl = 1:length(trialONTimes)
%     filePath = [boxPath filesep 'StoryDB' filesep 'p2fa_py3' filesep 'Recall_Task' filesep paths.patientPath filesep paths.sessPath filesep ['FRTranscript_trial' num2str(trl) '_TS.txt']];
%     RecallData.FRFULLTS{trl, 1} = Utilities.readinTransTextGrid(filePath);
%     % Adjust the timestampes so they are absolute timestamps
%     RecallData.FRFULLTS_Adj{trl, 1} = RecallData.FRFULLTS{trl, 1};
%     RecallData.FRFULLTS_Adj{trl, 1}(:, 2:3) = cellfun(@(x) {x+audioStartTimes(trl)}, RecallData.FRFULLTS{trl}(:, 2:3));
%     audioEndTimes(trl, 1) = cell2mat(RecallData.FRFULLTS_Adj{trl, 1}(end, 3));
% end

%% Set up each session

% P71 Sess 1
if isequal(paths.sessPath, 'Recall_Session_1_20201121')
    
    taskStruct = load([paths.basePath filesep paths.taskPath filesep paths.patientPath...
        filesep paths.sessPath filesep 'P71CS_1_Sub_']);
    
    % make encoding order for the encoding/cued recall sections
    rows = size(cell2mat(taskStruct.blockOrders), 1);
    cols = size(cell2mat(taskStruct.blockOrders), 2);
    EncodingOrder =  reshape(cell2mat(taskStruct.blockOrders)',...
        [rows*cols 1]);
    
    % order items were imagined in - taken from behavior notes and checked
    % with transcript
    CROrder = [1; 6; 1; 6; 1; 6;...
        4; 3; 4; 3; 4; 3;...
        6; 3; 6; 3; 6; 3;...
        5; 8; 5; 8; 5; 8;...
        8; 2; 8; 2; 8; 2;...
        7; 2; 7; 2; 7; 2;...
        7; 1; 7; 1; 7; 1;...
        4; 5; 4; 5; 4; 5];
    
    
    % Now find the absolute recall times for each image
    stimuli = {'Cat', 'Bird', 'Fruit', 'Oil', 'Cactus', 'Drops', 'Magnet', 'Shovel'}; % stimuli to be recalled
    RecallData.recEventMarkers(:, 1) = stimuli;
    RecallData.recEventMarkers(:, 2) = {1; 5; 2; 2; 4; 1; 6; 4}; % trials they were present in the first time
    RecallData.recEventMarkers(:, 3) = {[9; 19; 292]; [8; 39; 99]; [54; 109; 205]; [17; 97; 194]; [5; 20; 140]; [5; 73; 301]; [23; 70; 162]; [11; 26; 161]}; % markers within the trial they were recalled
    
    RecallData.recEventMarkers(:, 4) = {7; 6; 3; 8; 8; 3; 7; 5}; % trials they were present in the second time
    RecallData.recEventMarkers(:, 5) = {[11; 57; 166]; [7; 37; 134]; [8; 26; 82]; [3; 82; 153]; [9; 93; 156]; [16; 36; 90]; [7; 18; 151]; [23; 34; 87]}; % markers within the trial they were recalled
    
%     for stim = 1:length(stimuli)
%         RecallData.recEventMarkers{stim, 6} = [RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 2}}(RecallData.recEventMarkers{stim, 3}, 2);...
%             RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 4}}(RecallData.recEventMarkers{stim, 5}, 2)];
%         
%         % Check to make sure the right indices were chosen - these should
%         % be the spoken words of th relevant stimuli
%         RecallData.recEventMarkers{stim, 7} = [RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 2}}(RecallData.recEventMarkers{stim, 3}, 1);...
%             RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 4}}(RecallData.recEventMarkers{stim, 5}, 1)];
%     end
    
    
% P71 Sess 2
elseif isequal(paths.sessPath, 'Recall_Session_2_20201124')
    
    taskStruct = load([paths.basePath filesep paths.taskPath filesep paths.patientPath...
        filesep paths.sessPath filesep 'P71CS_Rec@_@_Sub_']);
    
    rows = size(cell2mat(taskStruct.blockOrders), 1);
    cols = size(cell2mat(taskStruct.blockOrders), 2);
    EncodingOrder =  reshape(cell2mat(taskStruct.blockOrders)',...
        [rows*cols 1]);
    
    % order items were imagined in - taken from behavior notes and checked
    % with transcript
    CROrder = [10; 5; 10; 5; 10; 5;...
        8; 9; 8; 9; 8; 9;...
        2; 5; 2; 5; 2; 5;...
        10; 9; 10; 9; 10; 9;...
        1; 12; 1; 12; 1; 12;...
        6; 7; 6; 7; 6; 7;...
        3; 6; 3; 6; 3; 6;...
        2; 3; 2; 3; 2; 3;...
        1; 4; 1; 4; 1; 4;...
        4; 7; 4; 7; 4; 7;...
        12; 11; 12; 11; 12; 11;...
        11; 8; 11; 8; 11; 8];
    
    stimuli = {'Head1', 'Jet', 'Head2', 'Truck', 'Head3', 'Biplane', 'Head4', 'Grapes', 'Head5', 'Vase', 'Head6', 'Lamp'}; % stimuli to be recalled
    
    RecallData.recEventMarkers(:, 1) = stimuli;
    RecallData.recEventMarkers(:, 2) = {5; 3; 7; 9; 1; 6; 6; 2; 2; 1; 11; 5}; % trials they were present in the first time
    RecallData.recEventMarkers(:, 3) = {[8; 54; 172]; [8; 24; 106]; [23; 83; 148]; [7; 21; 135]; [10; 38; 82]; [4; 26; 137];...
        [7; 79; 157]; [16; 27; 120]; [20; 67; 124]; [7; 14; 79]; [5; 56; 114]; [4; 14; 134]}; % markers within the trial they were recalled
    
    RecallData.recEventMarkers(:, 4) = {9; 8; 8; 10; 3; 7; 10; 12; 4; 4; 12; 11}; % trials they were present in the second time
    RecallData.recEventMarkers(:, 5) = {[13; 60; 143]; [18; 39; 146]; [23; 69; 165]; [2; 13; 85]; [10; 45; 123]; [18; 59; 143];...
        [6; 30; 94]; [8; 24; 99]; [11; 86; 165]; [6; 48; 139]; [12; 34; 105]; [2; 14; 100]}; % markers within the trial they were recalled
    
%     for stim = 1:length(stimuli)
%         RecallData.recEventMarkers{stim, 6} = [RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 2}}(RecallData.recEventMarkers{stim, 3}, 2);...
%             RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 4}}(RecallData.recEventMarkers{stim, 5}, 2)];
%         
%         % Check to make sure the right indices were chosen - these should
%         % be the spoken words of th relevant stimuli
%         RecallData.recEventMarkers{stim, 7} = [RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 2}}(RecallData.recEventMarkers{stim, 3}, 1);...
%             RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 4}}(RecallData.recEventMarkers{stim, 5}, 1)];
%     end

% P76 Sess 1
elseif isequal(paths.sessPath, 'ReScreenRecall_Session_1_20210917')
    
    taskStruct = load([paths.basePath filesep paths.taskPath filesep paths.patientPath...
        filesep paths.sessPath filesep 'P76CS_Rec_Sub_']);
    
    % make encoding order for the encoding/cued recall sections
    rows = size(cell2mat(taskStruct.blockOrders), 1);
    cols = size(cell2mat(taskStruct.blockOrders), 2);
    EncodingOrder =  reshape(cell2mat(taskStruct.blockOrders)',...
        [rows*cols 1]);
    
    % order items were imagined in - taken from behavior notes and checked
    % with transcript
    CROrder = [4; 5; 4; 5; 4; 5; 4; 5;...
        5; 2; 5; 2; 5; 2; 5; 2; ...
        4; 3; 4; 3; 4; 3; 4; 3;...
        1; 6; 1; 6; 1; 6; 1; 6;
        2; 3; 2; 3; 2; 3; 2; 3;...
        1; 6; 1; 6; 1; 6; 1; 6];
    
    
    % Now find the absolute recall times for each image
    stimuli = {'Cow', 'Puppy', 'Duckings', 'Microwave', 'LetterE', 'Stroller'}; % stimuli to be recalled
    RecallData.recEventMarkers(:, 1) = stimuli;
    RecallData.recEventMarkers(:, 2) = {4; 2; 3; 1; 1; 4}; % trials they were present in the first time
    RecallData.recEventMarkers(:, 3) = {[54; 212; 490]; [54; 137; 375]; [45; 201; 414]; [18; 39; 305]; [25; 83; 314]; [57; 294; 544]}; % markers within the trial they were recalled

     
    RecallData.recEventMarkers(:, 4) = {6; 5; 5; 3; 2; 6}; % trials they were present in the second time
    RecallData.recEventMarkers(:, 5) = {[41; 73; 468]; [39; 95; 476]; [45; 239; 481]; [24; 149; 393]; [39; 65; 358]; [38; 247; 590]}; % markers within the trial they were recalled

     
%     for stim = 1:length(stimuli)
%         RecallData.recEventMarkers{stim, 6} = [RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 2}}(RecallData.recEventMarkers{stim, 3}, 2);...
%             RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 4}}(RecallData.recEventMarkers{stim, 5}, 2)];
%         
%         % Check to make sure the right indices were chosen - these should
%         % be the spoken words of th relevant stimuli
%         RecallData.recEventMarkers{stim, 7} = [RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 2}}(RecallData.recEventMarkers{stim, 3}, 1);...
%             RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 4}}(RecallData.recEventMarkers{stim, 5}, 1)];
%     end

% P76 Sess 2
elseif isequal(paths.sessPath, 'Recall_Session_2_20210925')
    
    taskStruct = load([paths.basePath filesep paths.taskPath filesep paths.patientPath...
        filesep paths.sessPath filesep 'P76CSReScreen_Rec_2_Sub_']);
    
    % make encoding order for the encoding/cued recall sections
    rows = size(cell2mat(taskStruct.blockOrders), 1);
    cols = size(cell2mat(taskStruct.blockOrders), 2);
    EncodingOrder =  reshape(cell2mat(taskStruct.blockOrders)',...
        [rows*cols 1]);
    
    % order items were imagined in - taken from behavior notes and checked
    % with transcript
    % why is this identical to the previous days...
    CROrder = [4; 5; 4; 5; 4; 5; 4; 5; ...
        5; 2; 5; 2; 5; 2; 5; 2; ...
        4; 3; 4; 3; 4; 3; 4; 3; ...
        1; 6; 1; 6; 1; 6; 1; 6; ...
        3; 2; 3; 2; 3; 2; 3; 2; ...
        1; 6; 1; 6; 1; 6; 1; 6; ];
    
    
    % Now find the absolute recall times for each image
    stimuli = {'Parrot', 'HandWiphone', 'Tablets', 'Boy', 'LetterE', 'Pot'}; % stimuli to be recalled
    
    RecallData.recEventMarkers(:, 1) = stimuli;
    RecallData.recEventMarkers(:, 2) = {4; 2; 3; 1; 1; 4}; % trials they were present in the first time
    RecallData.recEventMarkers(:, 3) = {[41; 110; 393]; [37; 92; 239]; [26; 87; 211]; [68; 216; 411]; [127; 283; 491]; [92; 197; 333]}; % markers within the trial they were recalled
     
    RecallData.recEventMarkers(:, 4) = {6; 5; 5; 3; 2; 6}; % trials they were present in the second time
    RecallData.recEventMarkers(:, 5) = {[40; 46; 283]; [28; 147; 296]; [20; 94; 323]; [8; 36; 247]; [16; 44; 298]; [35; 141; 327]}; % markers within the trial they were recalled
     
%     for stim = 1:length(stimuli)
%         RecallData.recEventMarkers{stim, 6} = [RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 2}}(RecallData.recEventMarkers{stim, 3}, 2);...
%             RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 4}}(RecallData.recEventMarkers{stim, 5}, 2)];
%         
%         % Check to make sure the right indices were chosen - these should
%         % be the spoken words of th relevant stimuli
%         RecallData.recEventMarkers{stim, 7} = [RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 2}}(RecallData.recEventMarkers{stim, 3}, 1);...
%             RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 4}}(RecallData.recEventMarkers{stim, 5}, 1)];
%     end

% P76 Sess 3
elseif  isequal(paths.sessPath, 'ReScreenRecall_Session_3_20210927')
    
    taskStruct = load([paths.basePath filesep paths.taskPath filesep paths.patientPath...
        filesep paths.sessPath filesep 'P76CS_Recall_3_Sub_']);
    
    % make encoding order for the encoding/cued recall sections
    rows = size(cell2mat(taskStruct.blockOrders), 1);
    cols = size(cell2mat(taskStruct.blockOrders), 2);
    EncodingOrder =  reshape(cell2mat(taskStruct.blockOrders)',...
        [rows*cols 1]);
    
    % order items were imagined in - taken from behavior notes and checked
    % with transcript
    % why is this identical to the previous days...
    CROrder = [8; 3; 8; 3; 8; 3; 8; 3; ...
        5; 6; 5; 6; 5; 6; 5; 6; ...
        4; 5; 4; 5; 4; 5; 4; 5; ...
        6; 7; 6; 7; 6; 7; 6; 7; ...
        1; 2; 1; 2; 1; 2; 1; 2; ...
        1; 2; 1; 2; 1; 2; 1; 2; ...
        4; 3; 4; 3; 4; 3; 4; 3; ...
        8; 7; 8; 7; 8; 7; 8; 7; ];
    
    
    % Now find the absolute recall times for each image
    stimuli = {'Puppy', 'Lion', 'Lizard', 'Whale', 'Man', 'Woman', 'Pears', 'Drumsticks'}; % stimuli to be recalled
    
    RecallData.recEventMarkers(:, 1) = stimuli;
    RecallData.recEventMarkers(:, 2) = {5; 5; 1; 3; 2; 2; 4; 1}; % trials they were present in the first time
    RecallData.recEventMarkers(:, 3) = {[26; 30; 244]; [22; 70; 249]; [67; 206; 437]; [32; 131; 302];...
     [13; 25; 252]; [20; 83; 229]; [50; 181; 333]; [54; 86; 392]}; % markers within the trial they were recalled
    
    RecallData.recEventMarkers(:, 4) = {6; 6; 7; 7; 3; 4; 8; 8}; % trials they were present in the second time
    RecallData.recEventMarkers(:, 5) = {[130; 157; 318]; [133; 208; 374]; [10; 14; 281]; [6; 110; 230];...
    [14; 37; 324]; [31; 69; 374]; [22; 124; 256]; [8; 25; 302]}; % markers within the trial they were recalled
    
%     for stim = 1:length(stimuli)
%         RecallData.recEventMarkers{stim, 6} = [RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 2}}(RecallData.recEventMarkers{stim, 3}, 2);...
%             RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 4}}(RecallData.recEventMarkers{stim, 5}, 2)];
%         
%         % Check to make sure the right indices were chosen - these should
%         % be the spoken words of th relevant stimuli
%         RecallData.recEventMarkers{stim, 7} = [RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 2}}(RecallData.recEventMarkers{stim, 3}, 1);...
%             RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 4}}(RecallData.recEventMarkers{stim, 5}, 1)];
%     end

% P79 Sess 1
elseif isequal(paths.sessPath, 'ReScreenRecall_Session_1_20220330')
    
    % attempt 1  - Cedrus fucked me
    % Order reconstructed from audio file
%     CROrder = [5; 7; 5; 7; 5; 7; 5; 7;...
%           6; 3; 6; 3; 6; 3; 6; 3;...
%           7; 3; 7; 3; 7; 3; 7; 3;...
%           4; 8; 4; 8; 4; 8; 4; 8;...
%           4; 1; 4; 1; 4; 1; 4; 1;... % cedrus got unplugged here
%           ]
    
    taskStruct = load([paths.basePath filesep paths.taskPath filesep paths.patientPath...
        filesep paths.sessPath filesep 'P79CS_ReScreenRecall_1_Att2_Sub_']);
    
    rows = size(cell2mat(taskStruct.blockOrders), 1);
    cols = size(cell2mat(taskStruct.blockOrders), 2);
    EncodingOrder =  reshape(cell2mat(taskStruct.blockOrders)',...
        [rows*cols 1]);
    
    CROrder = [5; 6; 5; 6; 5; 6; 5; 6; ...
        7; 4; 7; 4; 7; 4; 7; 4; ...
        3; 7; 3; 7; 3; 7; 3; 7; ...
        4; 8; 4; 8; 4; 8; 4; 8; ...
        1; 2; 1; 2; 1; 2; 1; 2; ...
        1; 2; 1; 2; 1; 2; 1; 2; ...
        3; 6; 3; 6; 3; 6; 3; 6; ...
        5; 8; 5; 8; 5; 8; 5; 8; ];
     % Now find the absolute recall times for each image
    stimuli = {'Leopard', 'Brunette', 'Man looking right', 'Blonde', 'Man looking left', 'Shelves', 'Candle', 'Scope'}; % stimuli to be recalled
    
%     RecallData.recEventMarkers(:, 1) = stimuli;
%     RecallData.recEventMarkers(:, 2) = {5; 5; 1; 3; 2; 2; 4; 1}; % trials they were present in the first time
%     RecallData.recEventMarkers(:, 3) = {[26; 30; 244]; [22; 70; 249]; [67; 206; 437]; [32; 131; 302];...
%      [13; 25; 252]; [20; 83; 229]; [50; 181; 333]; [54; 86; 392]}; % markers within the trial they were recalled
%     
%     RecallData.recEventMarkers(:, 4) = {6; 6; 7; 7; 3; 4; 8; 8}; % trials they were present in the second time
%     RecallData.recEventMarkers(:, 5) = {[130; 157; 318]; [133; 208; 374]; [10; 14; 281]; [6; 110; 230];...
%     [14; 37; 324]; [31; 69; 374]; [22; 124; 256]; [8; 25; 302]}; % markers within the trial they were recalled
%     
%     for stim = 1:length(stimuli)
%         RecallData.recEventMarkers{stim, 6} = [RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 2}}(RecallData.recEventMarkers{stim, 3}, 2);...
%             RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 4}}(RecallData.recEventMarkers{stim, 5}, 2)];
%         
%         % Check to make sure the right indices were chosen - these should
%         % be the spoken words of th relevant stimuli
%         RecallData.recEventMarkers{stim, 7} = [RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 2}}(RecallData.recEventMarkers{stim, 3}, 1);...
%             RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 4}}(RecallData.recEventMarkers{stim, 5}, 1)];
%     end

% P79 Sess 2
elseif isequal(paths.sessPath, 'ReScreenRecall_Session_2_20220403')
    
    taskStruct = load([paths.basePath filesep paths.taskPath filesep paths.patientPath...
        filesep paths.sessPath filesep 'P79CS_ReScreenRecall_3_Sub_']);
    
    rows = size(cell2mat(taskStruct.blockOrders), 1);
    cols = size(cell2mat(taskStruct.blockOrders), 2);
    EncodingOrder =  reshape(cell2mat(taskStruct.blockOrders)',...
        [rows*cols 1]);
    
     CROrder = [4; 7; 4; 7; 4; 7; 4; 7; ...
        5; 6; 5; 6; 5; 6; 5; 6; ...
        2; 6; 2; 6; 2; 6; 2; 6; ...
        5; 1; 5; 1; 5; 1; 5; 1; ...
        3; 7; 3; 7; 3; 7; 3; 7; ...
        8; 4; 8; 4; 8; 4; 8; 4; ...
        8; 3; 8; 3; 8; 3; 8; 3; ...
        1; 2; 1; 2; 1; 2; 1; 2; ];
    
    stimuli = {'Leopard', 'Cow', 'Flashlight', 'Toilet paper', 'Rainbow', 'Sun rays', 'Lightbulb', 'Binoculars'}; % stimuli to be recalled
    
    %     RecallData.recEventMarkers(:, 1) = stimuli;
%     RecallData.recEventMarkers(:, 2) = {5; 5; 1; 3; 2; 2; 4; 1}; % trials they were present in the first time
%     RecallData.recEventMarkers(:, 3) = {[26; 30; 244]; [22; 70; 249]; [67; 206; 437]; [32; 131; 302];...
%      [13; 25; 252]; [20; 83; 229]; [50; 181; 333]; [54; 86; 392]}; % markers within the trial they were recalled
%     
%     RecallData.recEventMarkers(:, 4) = {6; 6; 7; 7; 3; 4; 8; 8}; % trials they were present in the second time
%     RecallData.recEventMarkers(:, 5) = {[130; 157; 318]; [133; 208; 374]; [10; 14; 281]; [6; 110; 230];...
%     [14; 37; 324]; [31; 69; 374]; [22; 124; 256]; [8; 25; 302]}; % markers within the trial they were recalled
%     
%     for stim = 1:length(stimuli)
%         RecallData.recEventMarkers{stim, 6} = [RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 2}}(RecallData.recEventMarkers{stim, 3}, 2);...
%             RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 4}}(RecallData.recEventMarkers{stim, 5}, 2)];
%         
%         % Check to make sure the right indices were chosen - these should
%         % be the spoken words of th relevant stimuli
%         RecallData.recEventMarkers{stim, 7} = [RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 2}}(RecallData.recEventMarkers{stim, 3}, 1);...
%             RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 4}}(RecallData.recEventMarkers{stim, 5}, 1)];
%     end

% P79 Sess 3
elseif isequal(paths.sessPath, 'ReScreenRecall_Session_3_20220405')
    
    taskStruct = load([paths.basePath filesep paths.taskPath filesep paths.patientPath...
        filesep paths.sessPath filesep 'P79CS_ReScreenRecall_4_Sub_']);
    
    rows = size(cell2mat(taskStruct.blockOrders), 1);
    cols = size(cell2mat(taskStruct.blockOrders), 2);
    EncodingOrder =  reshape(cell2mat(taskStruct.blockOrders)',...
        [rows*cols 1]);
    
     CROrder = [7; 2; 7; 2; 7; 2; 7; 2; ...
        6; 3; 6; 3; 6; 3; 6; 3; ...
        1; 8; 1; 8; 1; 8; 1; 8; ...
        4; 1; 4; 1; 4; 1; 4; 1; ...
        8; 7; 8; 7; 8; 7; 8; 7; ...
        2; 5; 2; 5; 2; 5; 2; 5; ...
        4; 5; 4; 5; 4; 5; 4; 5; ...
        6; 3; 6; 3; 6; 3; 6; 3; ];
    
    stimuli = {'Turtle', 'Battery', 'Man', 'Blackberries', 'P', 'Piano', 'Palm Tree', 'Ribbon'}; % stimuli to be recalled
    
    %     RecallData.recEventMarkers(:, 1) = stimuli;
%     RecallData.recEventMarkers(:, 2) = {5; 5; 1; 3; 2; 2; 4; 1}; % trials they were present in the first time
%     RecallData.recEventMarkers(:, 3) = {[26; 30; 244]; [22; 70; 249]; [67; 206; 437]; [32; 131; 302];...
%      [13; 25; 252]; [20; 83; 229]; [50; 181; 333]; [54; 86; 392]}; % markers within the trial they were recalled
%     
%     RecallData.recEventMarkers(:, 4) = {6; 6; 7; 7; 3; 4; 8; 8}; % trials they were present in the second time
%     RecallData.recEventMarkers(:, 5) = {[130; 157; 318]; [133; 208; 374]; [10; 14; 281]; [6; 110; 230];...
%     [14; 37; 324]; [31; 69; 374]; [22; 124; 256]; [8; 25; 302]}; % markers within the trial they were recalled
%     
%     for stim = 1:length(stimuli)
%         RecallData.recEventMarkers{stim, 6} = [RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 2}}(RecallData.recEventMarkers{stim, 3}, 2);...
%             RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 4}}(RecallData.recEventMarkers{stim, 5}, 2)];
%         
%         % Check to make sure the right indices were chosen - these should
%         % be the spoken words of th relevant stimuli
%         RecallData.recEventMarkers{stim, 7} = [RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 2}}(RecallData.recEventMarkers{stim, 3}, 1);...
%             RecallData.FRFULLTS_Adj{RecallData.recEventMarkers{stim, 4}}(RecallData.recEventMarkers{stim, 5}, 1)];
%     end  

% P80 Sess 1
elseif isequal(paths.sessPath, 'ReScreenRecall_Session_1_20220728')
    
    taskStruct = load([paths.basePath filesep paths.taskPath filesep paths.patientPath...
        filesep paths.sessPath filesep 'P80CS_ReScreenRecall_Sub_']);
    
    rows = size(cell2mat(taskStruct.blockOrders), 1);
    cols = size(cell2mat(taskStruct.blockOrders), 2);
    EncodingOrder =  reshape(cell2mat(taskStruct.blockOrders)',...
        [rows*cols 1]);
    
     CROrder = [2; 6; 2; 6; 2; 6; 2; 6; ...
        3; 7; 3; 7; 3; 7; 3; 7; ...
        1; 5; 1; 5; 1; 5; 1; 5; ...
        7; 4; 7; 4; 7; 4; 7; 4; ...
        8; 3; 8; 3; 8; 3; 8; 3; ...
        6; 4; 6; 4; 6; 4; 6; 4; ...
        1; 5; 1; 5; 1; 5; 1; 5; ...
        8; 2; 8; 2; 8; 2; 8; 2; ];
    
    stimuli = {'Deer', 'Penguin', 'Tigers', 'Coffeemaker', 'Smiling lady', 'Asian guy', 'Lady', 'Boat'}; % stimuli to be recalled
    
    
% P80 Sess 2
% need to figure out - what stimuli I used for recall ( can I do this with marked positions?)
% CR block order 
elseif isequal(paths.sessPath, 'ReScreenRecall_Session_2_20220731')
    
    taskStruct = load([paths.basePath filesep paths.taskPath filesep paths.patientPath...
        filesep paths.sessPath filesep 'P80CS_Recall_2_Sub_']);
    
    rows = size(cell2mat(taskStruct.blockOrders), 1);
    cols = size(cell2mat(taskStruct.blockOrders), 2);
    EncodingOrder =  reshape(cell2mat(taskStruct.blockOrders)',...
        [rows*cols 1]);
    
     CROrder = [2; 6; 2; 6; 2; 6; 2; 6; ...
        3; 7; 3; 7; 3; 7; 3; 7; ...
        5; 1; 5; 1; 5; 1; 5; 1; ...
        4; 7; 4; 7; 4; 7; 4; 7; ...
        3; 8; 3; 8; 3; 8; 3; 8; ...
        4; 6; 4; 6; 4; 6; 4; 6; ...
        1; 5; 1; 5; 1; 5; 1; 5; ...
        8; 2; 8; 2; 8; 2; 8; 2; ];
    
    stimuli = {'Peacock', 'Cap', 'Man', 'Umbrella', 'Bee', 'L', 'Clown with drum', 'Log'}; % stimuli to be recalled
end
RecallData.stimuli = stimuli;
RecallData.EncodingOrder = EncodingOrder;
RecallData.CROrder = CROrder;
RecallData.trialONTimes = trialONTimes;
RecallData.trialOFFTimes = trialOFFTimes;
RecallData.recallONTimes = recallONTimes;
RecallData.recallOFFTimes = recallOFFTimes;
RecallData.toneONTimes = toneONTimes;
RecallData.toneOFFTimes = toneOFFTimes;
RecallData.imageONTimes = imageONTimes;
RecallData.imageOFFTimes = imageOFFTimes;
RecallData.audioStartTimes = audioStartTimes;
RecallData.audioEndTimes = audioEndTimes;
RecallData.distONTimes = distractionONTimes;
RecallData.distOFFTimes = distractionOFFTimes;
% RecallData.FRTimes = FRTimes;
end