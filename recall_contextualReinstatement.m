dbstop if error
[~, host] = system('hostname');
if strcmp(host(1:end-1), 'DWA644201')
    atCedars = 1;
    diskPath = 'G:\SUAnalysis';
elseif strcmp(host(1:end-1), 'DESKTOP-LJHLIED')
    atCedars = 0;
    diskPath = 'G:\SUAnalysis';
elseif strcmp(host(1:end-1), 'Varuns-iMac-2.local')
    atCedars = 0;
    diskPath = '/Volumes/T7/SUAnalysis';
end
taskPath = [diskPath filesep 'Recall_Task'];

Sess = 1;

if Sess == 1
    basePath = 'G:\SUAnalysis\Recall_Task\P76CS\ReScreenRecall_Session_1_20210917';    
elseif Sess == 2
    basePath = 'G:\SUAnalysis\Recall_Task\P76CS\Recall_Session_2_20210925'; 
elseif Sess == 3
    basePath = 'G:\SUAnalysis\Recall_Task\P76CS\ReScreenRecall_Session_3_20210927';
end


pathFrames = [basePath filesep 'stimuliUsedRecall'];


if strcmp(host(1:end-1), 'DWA644201')
    taskCodePath = 'D:\Users\wadiav\Dropbox\Caltech\Thesis\Human_work\Cedars\RecallTaskVarun';
elseif strcmp(host(1:end-1), 'DESKTOP-LJHLIED')
    taskCodePath = 'E:\Dropbox\Caltech\Thesis\Human_work\Cedars\RecallTaskVarun';
elseif strcmp(host(1:end-1), 'Varuns-iMac-2.local')
    taskCodePath = '/Users/varunwadia/Dropbox/Caltech/Thesis/Human_work/Cedars/RecallTaskVarun';
end


addpath(taskCodePath);

setTTLCodes


%% load in appropriate RecallData struct

AllCells = 0; IT_MTL = 0; ITOnly = 1;

if AllCells
    load([basePath filesep 'RecallData']);
    load([basePath filesep 'PsthandResponses']);
elseif IT_MTL
    load([basePath filesep 'RecallData_IT_MTL']);
    load([basePath filesep 'IT_MTL_PsthandResponses.mat']);
elseif ITOnly
    load([basePath filesep 'RecallData_IT']);
    load([basePath filesep 'ITPsthandResponses']);
end


%%
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




%% creating encoding vec and CR matrix (no-zscoring) per stim

% Need to create population vector for each stim during screening/encoding (depending on the region)
% Split the CR sections into bins and correlate the pop vec with each time step and plot that correlation value 
% Do this both per trial (Like Jie did - so your output *per stim* is a matrix cells*trials x bins) and with the trial averages (so matrix per stim is cells x bins)
% Also make line plot for each case.
% Do I zscore each cell's response??


% compute encoding spike count using respLat from screeing OR fixed time
% bin if that doesn't exist

normalize = 0;

offset = 100; %ms
windowLength = RecallData.stimDur*1e3;

[sortedOrder, correctOrder] = sortrows(RecallData.EncodingOrder);
EncResp_perTrial = [];

stimCorrMat = zeros(length(RecallData.stimuli), 1);
for stim = 1:length(RecallData.stimuli)
    if normalize
        EncMat = zscore(RecallData.perStimAllCellsEncoding{stim, 1}')';
        EncResp_perTrial = mean(EncMat(:, (-RecallData.timelimits(1)*1e3)+offset:ceil((-RecallData.timelimits(1)*1e3)+offset+windowLength)), 2);
        
        CRResp_perTrial = zscore(RecallData.perStimAllCellsCR{stim, 1}')';
        
    else
        % if I use fixed time window for all cells
        EncResp_perTrial = mean(RecallData.perStimAllCellsEncoding{stim, 1}(:, (-RecallData.timelimits(1)*1e3)+offset:ceil((-RecallData.timelimits(1)*1e3)+offset+windowLength)), 2); % response per trial
        CRResp_perTrial = RecallData.perStimAllCellsCR{stim, 1}; % this is already sorted
    end
    ctr = 1;
    for bin = 1:100:length(RecallData.CRTimeCourse{1, 3})-floor(windowLength)
        corr_specificTimeBin = corrcoef(EncResp_perTrial, mean(CRResp_perTrial(:, bin:bin+floor(windowLength)), 2));
        stimCorrMat(stim, ctr) = corr_specificTimeBin(1, 2);
        ctr = ctr+1;
    end
end

%% Shuffle distribution

% find the time point for which contREin is peaked 
% for that time bin shuffle the labels and form distribution - see where it
% lies

[peaks, t_step] = max(stimCorrMat(:, 10:end), [], 2);
t_step = t_step+9;

n_reps = 1000;
p_vals = ones(length(t_step), 1);

for t = 1:length(t_step)
    
    ShuffCorr = [];

    
    for n = 1:n_reps
        
        if t_step(t)*100+floor(windowLength) > 7500
            cc = corrcoef(EncResp_perTrial, Shuffle(mean(CRResp_perTrial(:, t_step(t)*100:end), 2))); 
        else
            cc = corrcoef(EncResp_perTrial, Shuffle(mean(CRResp_perTrial(:, t_step(t)*100:t_step(t)*100+floor(windowLength)), 2))); 
        end
        ShuffCorr = [ShuffCorr cc(1,2)];
        
    end
    p_vals(t) = sum(ShuffCorr > peaks(t))/length(ShuffCorr);
end


%% plot


f = figure; 
% hold on
imagesc(stimCorrMat);
% imagesc(stimCorrMat, [0 max(max(stimCorrMat))]);
% setting up my colormap
orangemap = esa(300);
[WhiteColor, Whitepos] = max(sum(orangemap,2));
orangemap = orangemap([1:Whitepos round(linspace(Whitepos+1,size(orangemap,1)-2,Whitepos-1))],:);
colormap(orangemap)
colbar = colorbar;
colbar.AxisLocation = 'out';
colbar.Position = [0.92,0.70,0.020,0.22];
% colbar.Label.String = 'Corr coef';

if normalize
    if AllCells
        title('Contextual reinstatement during imagination period')
        filename = [basePath filesep 'ContextRein_AllCells_Norm_Sess_' num2str(Sess) '_forPaper'];
    elseif IT_MTL
        title('Contextual reinstatement of IT+MTL population during imagination period')
        filename = [basePath filesep 'ContextRein_IT_MTL_Norm_Sess_' num2str(Sess) '_forPaper'];
        
    elseif ITOnly
        title('Contextual reinstatement of IT population during imagination period')
        filename = [basePath filesep 'ContextRein_IT_Norm_Sess_' num2str(Sess) '_forPaper'];     
    end
    
else
    if AllCells
        title('Contextual reinstatement during imagination period')
        filename = [basePath filesep 'ContextRein_AllCells_Sess_' num2str(Sess) '_forPaper'];
    elseif IT_MTL
        title('Contextual reinstatement of IT+MTL population during imagination period')
        filename = [basePath filesep 'ContextRein_IT_MTL_Sess_' num2str(Sess) '_forPaper'];
        
    elseif ITOnly
        title('Contextual reinstatement of IT population during imagination period')
        filename = [basePath filesep 'ContextRein_IT_Sess_' num2str(Sess) '_forPaper'];
    end
end

xlabel('Time (s)');
xticks([1 11 60]);

set(gca, 'xticklabel', {[]}); % clear the old labels
xtik = get(gca,'xtick'); 
xticklabels({'-0.5', '0', '5'});

yl = ylim;
hold on
plot([11 11], [yl(1) yl(2)], '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
set(gca, 'FontWeight', 'bold');

print(f, filename, '-dpng', '-r0')

% if I wanted to use a different respLat for each cell - not necessary most likely
% for cellIndex = l(strctCells)
%     
%     if ~isempty(responses{cellIndex, 2})
%         respLat = (-RecallData.timelimits(1)*1e3)+responses{cellIndex, 2};
%     else
%         respLat = (-RecallData.timelimits(1)*1e3)+offset+RecallData.stimDur*1e3;
%     end
%     
%     
%     CellEncResp = mean(RecallData.EncodingTimeCourse{cellIndex, 1}(:, respLat:floor(respLat+windowLength)), 2);
%     EncResp_perTrial = [EncResp_perTrial; CellEncResp(correctOrder, :)];
%     
% end

%% xcorr

% take two vectors and feed them into xcorr function
% According to Hristos: 
% xcorr after smoothing with a kernel is one way
% computing a cross-correlogram directly is another way
% (ccg = go through all spikes in one spike train and compute the peri-spike firing rate histogram using the other spike train)


% Take smoothed psth for a cell in IT and a cell in Hipp and compute.
% do this for all combinations of cells 















