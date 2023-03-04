% trialNum: N x S x D
% firingRates: N x S x D x T x maxTrialNum
% firingRatesAverage: N x S x D x T
%
% N is the number of neurons
% S is the number of stimuli conditions (boundary types in Jie's task/ encoding, CR and FR in my task)
% D is the number of decisions (D=2) % not relavant in my task
% --------------------------------------------------------------------
% T is the number of time-points (note that all the trials should have the
% same length in time!)
% --------------------------------------------------------------------
% trialNum -- number of trials for each neuron in each S,D condition (is
% usually different for different conditions and different sessions)
%
% firingRates -- all single-trial data together, massive array. Here
% maxTrialNum is the maximum value in trialNum. E.g. if the number of
% trials per condition varied between 1 and 20, then maxTrialNum = 20. For
% the neurons and conditions with less trials, fill remaining entries in
% firingRates with zeros or nans.
%
% firingRatesAverage -- average of firingRates over trials (5th dimension).
% If the firingRates is filled up with nans, then it's simply
%    firingRatesAverage = nanmean(firingRates,5)
% If it's filled up with zeros (as is convenient if it's stored on hard 
% drive as a sparse matrix), then 
%    firingRatesAverage = bsxfun(@times, mean(firingRates,5), size(firingRates,5)./trialNum)

%% compressing population responses - for stateSpace analyses
% Need a neurons x time matrix for both CR and FR
% trialNum: n_cells x n_CondTypes (1: CR; 2: FR)
% firingRates: n_cells x n_CondTypes x n_times x n_maxTrial
% firingRatesAverage: nanmean(firintRates, 4), n_cells x n_BTypes x n_times

timePeriodCR = length(RecallData.perStimAllCellsCR{1, 3});
timePeriodFR = length(RecallData.perStimAllCellsFR{1, 3});
assert(isequal(timePeriodCR, timePeriodFR)); % these need to have the same time

n_cells = length(strctCells);
n_CondTypes = 2; % 1 CR and 2 FR
n_maxTrial = length(RecallData.trialONTimes); % maximal number of trial repetitions


trialNum = nan(n_cells, n_CondTypes);

withoutBinning = 0;
binWidth = 10;
if withoutBinning
    n_times = timePeriodCR; % number of time points
else
    n_times = timePeriodCR/binWidth;
end
trialRaster_CR = nan(length(strctCells), timePeriodCR);
trialRaster_CR_norm = nan(length(strctCells), n_times);
trialRaster_FR = nan(length(strctCells), timePeriodCR);
trialRaster_FR_norm = nan(length(strctCells), n_times);
firingRates = nan(n_cells, n_CondTypes, n_times, n_maxTrial);


for trl = 1:length(RecallData.trialONTimes)
    for cellIndex = 1:length(strctCells)
        % without binning take the smoothed psth directly
        if withoutBinning
            n_times = timePeriodCR;
            trialRaster_CR(cellIndex, :) = mean(RecallData.CRTimeCourse{cellIndex, 2}(find(RecallData.CROrder == trl), :), 1);
            trialRaster_FR(cellIndex, :) = mean(RecallData.perStimAllCellsFR{trl, 2}(find(RecallData.order_perStimAllCellsFR == cellIndex), :), 1);
            tR_CR(cellIndex, :) = trialRaster_CR(cellIndex, :);
            tR_FR(cellIndex, :) = trialRaster_FR(cellIndex, :);
        else
            n_times = timePeriodCR/binWidth;

            % if need to bin first take raster
            trialRaster_CR(cellIndex, :) = mean(RecallData.CRTimeCourse{cellIndex, 1}(find(RecallData.CROrder == trl), :), 1);
            trialRaster_FR(cellIndex, :) = mean(RecallData.perStimAllCellsFR{trl, 1}(find(RecallData.order_perStimAllCellsFR == cellIndex), :), 1);
            ctr = 1;
            for bin = 1:binWidth:timePeriodCR
                % bin
                trialRaster_CR_binned(cellIndex, ctr) = sum(trialRaster_CR(cellIndex, bin:bin+binWidth-1)); 
                trialRaster_FR_binned(cellIndex, ctr) = sum(trialRaster_FR(cellIndex, bin:bin+binWidth-1)); 
                
                ctr = ctr + 1;
            end
            % smooth
            trialRaster_CR_binned(cellIndex, :) = smoothdata(trialRaster_CR_binned(cellIndex, :), 2, 'gaussian', 20);
            trialRaster_FR_binned(cellIndex, :) = smoothdata(trialRaster_FR_binned(cellIndex, :), 2, 'gaussian', 20);
            tR_CR(cellIndex, :) = trialRaster_CR_binned(cellIndex, :);
            tR_FR(cellIndex, :) = trialRaster_FR_binned(cellIndex, :);
            
            
        end
        
        
        % normalize
        trialMean_CR = mean(tR_CR(cellIndex, :), 2);
        trialStd_CR = std(tR_CR(cellIndex, :), [], 2);
        trialMean_FR = mean(tR_FR(cellIndex, :), 2);
        trialStd_FR = std(tR_FR(cellIndex, :), [], 2);
        
        if trialStd_CR == 0 && trialStd_FR ~= 0
            trialRaster_CR_norm(cellIndex, :) = zeros(1, n_times);
            trialRaster_FR_norm(cellIndex, :) = (tR_FR(cellIndex, :) - repmat(trialMean_FR, 1, n_times))./(repmat(trialStd_FR, 1, n_times));
        elseif trialStd_FR == 0 && trialStd_CR ~= 0
            trialRaster_CR_norm(cellIndex, :) = (tR_CR(cellIndex, :) - repmat(trialMean_CR, 1, n_times))./(repmat(trialStd_CR, 1, n_times));
            trialRaster_FR_norm(cellIndex, :) = zeros(1, n_times);
        elseif trialStd_FR == 0 && trialStd_CR == 0
            trialRaster_CR_norm(cellIndex, :) = zeros(1, n_times);
            trialRaster_FR_norm(cellIndex, :) = zeros(1, n_times);
        else
            % normalize
            trialRaster_CR_norm(cellIndex, :) = (tR_CR(cellIndex, :) - repmat(trialMean_CR, 1, n_times))./(repmat(trialStd_CR, 1, n_times));
            trialRaster_FR_norm(cellIndex, :) = (tR_FR(cellIndex, :) - repmat(trialMean_FR, 1, n_times))./(repmat(trialStd_FR, 1, n_times));
        end
        
        firingRates(cellIndex, 1, :, trl) = trialRaster_CR_norm(cellIndex, :);
        firingRates(cellIndex, 2, :, trl) = trialRaster_FR_norm(cellIndex, :);
        %         RecallData.CR_stateSpace(cellIndex, 1, :, trl) = mean(RecallData.perStimAllCellsCR{trl, 1}(find(RecallData.order_perStimAllCellsCR == cellIndex), :), 1);
        %         RecallData.FR_stateSpace(cellIndex, 1, :, trl) = mean(RecallData.perStimAllCellsFR{trl, 1}(find(RecallData.order_perStimAllCellsFR == cellIndex), :), 1);
        trialNum(cellIndex, 1) = 8;
        trialNum(cellIndex, 2) = 8;
    end
end

firingRatesAverage = nanmean(firingRates, 4);

%% Define parameter grouping

% For two parameters (e.g. stimulus and time, but no decision), we would have
% firingRates array of [N S T E] size (one dimension less, and only the following
% possible marginalizations:
%    1 - stimulus
%    2 - time
%    [1 2] - stimulus/time interaction
% They could be grouped as follows: 
%    combinedParams = {{1, [1 2]}, {2}};

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines

combinedParams = {{1, [1 2]}, {2}};
margNames = {'CondType', 'Time'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

time = RecallData.CRTimeCourse{1, 3}./1e3; % the length of a trial in seconds
timeEvents = 0; % the time a boundary happens in her case, time a recall event begins in mine 

% check consistency between trialNum and firingRates
for n = 1:size(firingRates,1)
    for s = 1:size(firingRates,2)
            assert(isempty(find(isnan(firingRates(n,s,:,1:trialNum(n,s))), 1)), 'Something is wrong!')
    end
end

%% PCA pf dataset

X = firingRatesAverage(:,:);
X = bsxfun(@minus, X, mean(X,2));

% -----------------------
[W,~,~] = svd(X, 'econ');

% computing explained variance
explVar = dpca_explainedVariance(firingRatesAverage, W, W, ...
    'combinedParams', combinedParams); 
% --------------------------
% tic
% [W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
%     'combinedParams', combinedParams);
% toc
% 
% explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
%     'combinedParams', combinedParams);
% --------------------------

% loadings in each PC
X = firingRatesAverage(:,:)';
Xcen = bsxfun(@minus, X, mean(X));
Z = Xcen * W;
componentsToPlot = find(explVar.cumulativePCA < 100);
dataDim = size(firingRatesAverage);
Zfull = reshape(Z(:,componentsToPlot)', [length(componentsToPlot) dataDim(2:end)]); % n_PCs x n_BTypes x n_times
CR_PC = squeeze(Zfull(:, 1, :));
FR_PC = squeeze(Zfull(:, 2, :)); % am I doing this right?


%% plot the top 3 PCs

CR_color = [255,99,71]./255;
FR_color = [100,149,237]./255;

figure('rend','painters','pos',[10 10 450 400]) % CR
if withoutBinning
    plot3(CR_PC(1,100:end-100),CR_PC(2,100:end-100),CR_PC(3,100:end-100),  'LineWidth', 4, 'color', CR_color);
else
    plot3(CR_PC(1,:),CR_PC(2,:),CR_PC(3,:),  'LineWidth', 4, 'color', CR_color);
end
hold on % if placed earlier it freezes the axes properties and makes 2d plot
plot3(timeEvents, timeEvents, timeEvents, '.', 'MarkerSize', 30, 'color', 'black')
xlim([-2 2]); xlabel('PC1');
ylim([-2 2]); ylabel('PC2');
zlim([-2 2]); zlabel('PC3');
title('Cued Recall');
set(gca, 'LineWidth', 1.5, 'FontSize', 20, 'FontWeight', 'bold');
box on
grid on

figure('rend','painters','pos',[10 10 450 400]) % FR
if withoutBinning
    plot3(FR_PC(1,100:end-100),FR_PC(2,100:end-100),FR_PC(3,100:end-100),  'LineWidth', 4, 'color', FR_color);
else
    plot3(FR_PC(1,:),FR_PC(2,:),FR_PC(3,:),  'LineWidth', 4, 'color', FR_color); 
end
hold on
plot3(timeEvents, timeEvents, timeEvents, '.', 'MarkerSize', 30, 'color', 'black')
xlim([-2 2]); xlabel('PC1');
ylim([-2 2]); ylabel('PC2');
zlim([-2 2]); zlabel('PC3');
title('Free Recall');
set(gca, 'LineWidth', 1.5, 'FontSize', 20, 'FontWeight', 'bold');
box on
grid on


