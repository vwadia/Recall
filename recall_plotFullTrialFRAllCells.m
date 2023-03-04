clear plotOptions
clear textParams

% strctCELL = struct2cell(strctCells');
% strctCELL = strctCELL';
cellNames = strctCELL(:, 1);

textParams.xlabel = 'time(MS)';
textParams.ylabelplot2 = 'trials (re-ordered)';
textParams.ylabelplot1 = 'Firing Rate (Hz)';
textParams.legendTitle = 'Cells';

popOrder = cell2mat(RecallData.order_perTrialAllCellsFullFR);
brainAreas = unique(popOrder, 'stable');
regionNames = unique(strctCELL(:, 4), 'stable');

isPDF = 0;


for trl = 1:length(RecallData.trialONTimes)
    
    FROrder = popOrder; % per region
    FRCols = Utilities.distinguishable_colors(length(brainAreas));
    
    textParams.fold_name = ['PopulationFRResponse'];
    textParams.suptitle = ['Population response to FRPeriod trial ' num2str(trl)];
    
    times = RecallData.perTrialAllCellsFullFR{trl, 3};
    transcription = RecallData.FRFULLTS{trl, 1};
    
    % legend
    areaCellNames = regionNames;
    for ar = 1:length(areaCellNames)
        areaCellNames{ar} = num2str(areaCellNames{ar});
    end
    textParams.legend = areaCellNames;
    %         textParams.useTheseColors = fullCols(((splt-1)*length(regCells))+1:(splt*length(regCells)), :);
    paths.destPathRecall = [paths.basePath filesep paths.taskPath filesep paths.patientPath filesep paths.sessPath filesep...
        'processedData' filesep textParams.fold_name];
    if ~exist(paths.destPathRecall)
        mkdir(paths.destPathRecall);
    end
    raster = RecallData.perTrialAllCellsFullFR{trl, 1};
    specificColor = [];
    handlesToFig = Utilities.Plotting.plotRastersTextinChunks_2(paths, raster, times, FROrder, RecallData.offsetFR, transcription, textParams, specificColor, isPDF);%, strctCells);
    
end