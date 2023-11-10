% script to concatenate all IT cells in recall (for plotting purposes)


setDiskPaths

taskCodePath = [boxPath filesep 'RecallTaskVarun'];
addpath(taskCodePath)
setTTLCodes

% set session list
task = 'Recall_Task';
sessID = Utilities.sessionListAllTasks(task, false);
sessID{2} = ['Recall_Task' filesep 'P76CS' filesep 'Recall_Session_2_20210925']; % want the rasters so using recall session

alpha = 0.05;

% load([diskPath filesep task filesep 'ReactivationCorrelationData0.05.mat'])
% RC_names = reacData.Name;
% RC_names = reacData.Channel;

% load([diskPath filesep task filesep ['ReactiveITCells_alpha' num2str(alpha) '_500Stim_Im_SigRamp.mat']]); 
% RC = strctCells;
% clearvars psths responses strctCells

% % load manually curated reactive cells
% RC_names_man = readtable([diskPath filesep task filesep 'ReactiveITCells_manual']);
% RC_nm = [table2array(RC_names_man(:, 1)), table2array(RC_names_man(:, 3))];
% RC_names = RC_nm(:, 1);
% RC_ChanNum = RC_nm(:, 2);

% RC_names = cat(1, RC(:).Name);
% RC_ChanNum = cat(1, RC(:).ChannelNumber);

% load([diskPath filesep task filesep  'im_active_cells_inThesis.mat'])
% RC_names = cat(1, cell2mat(im_active_cells(:, 1)));
% RC_ChanNum = cat(1, cell2mat(im_active_cells(:, 2)));

% RC_names = table2array(reacData(reacData.reactivated, 1));
% RC_ChanNum = table2array(reacData(reacData.reactivated, 2));

numReacTrialsPerStim = 3;
MinSpInBurst = 5;
load([diskPath filesep task filesep 'ReacDataPoisson_' num2str(numReacTrialsPerStim) 'ReacTrPerStim_MinSPB' num2str(MinSpInBurst)])
RC_names = reacData.Name(reacData.reactivated);
RC_ChanNum = reacData.Channel(reacData.reactivated);

% smoothingBin = 500; % for enc and CR - taken from recallscript_main

ITOnly = true;

plotOptions.MarkerSize = 4;
plotOptions.Fontsize = 16;
plotOptions.LineWidth = 3;
plotOptions.NormPsthFR = false;
plotOptions.fullscreen = true;
plotOptions.visible = false;
plotOptions.legend = true;
plotOptions.showReacTrials = true;
plotOptions.binSpikes = true;
plotOptions.binSize = 50; % ms
plotOptions.stepSize = 25;
plotOptions.smoothSize = 10;
% list of nicer colormaps
[mcs,nmn,pyt] = brewermap('list');


tic 
for ss = 1:length(sessID)
     
     
    load([sessID{ss} filesep 'RecallData_NoFreeRec.mat'])
    strctCELL = struct2cell(strctCells')';
    
    if ITOnly
        IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
        full_strctCELL = strctCELL;
        strctCELL = strctCELL(IT_Cells, :);
        strctCells = strctCells(IT_Cells);
        EncodingTimeCourse =  RecallData.EncodingTimeCourse(IT_Cells, 1:3);
        CRTimeCourse = RecallData.CRTimeCourse(IT_Cells, 1:3);
        
        % Baseline for CR
        BTimeCourse = RecallData.PreTrialBaselineTimeCourse(IT_Cells, 1:3);        
    else
        EncodingTimeCourse =  RecallData.EncodingTimeCourse(:, 1:3);
        CRTimeCourse = RecallData.CRTimeCourse(:, 1:3);
        
        % Baseline for CR
        BTimeCourse = RecallData.PreTrialBaselineTimeCourse(:, 1:3);                
    end
    
    num_cells(ss) = length(strctCells);

   
    for cellIndex = 1:length(strctCells)
        
        if ss == 1
            idx = cellIndex;
        else
            idx = sum(num_cells(1:ss-1)) + cellIndex;
        end
%     parfor cellIndex = 1:length(strctCells)
        enc_psth = EncodingTimeCourse(cellIndex, 1:3);
        CR_psth = CRTimeCourse(cellIndex, 1:3);
        
        if plotOptions.binSpikes
            
%             psthB = Utilities.Plotting.BinRasterSpikes(enc_psth{1, 1}, plotOptions.binSize/2, plotOptions.stepSize/2);
%             psthBS = [];
%             % smooth after binning
%             for row = 1:size(psthB, 1)
%                 psthBS(row, :) =  (Utilities.Smoothing.fastsmooth(psthB(row,:),plotOptions.smoothSize, 1, 0)*1e3);
%             end
%             enc_psth_binned{1, 2} = psthBS;
            
%             % smooth before binning
            enc_psth_binned{1, 2} = Utilities.Plotting.BinRasterSpikes(enc_psth{1, 2}, plotOptions.binSize, plotOptions.stepSize);
            
            if isfield(plotOptions, 'stepSize')
              enc_psth_binned{1, 3} = enc_psth{1, 3}(1:plotOptions.stepSize:end);  
            else
                enc_psth_binned{1, 3} = enc_psth{1, 3}(1:plotOptions.binSize:end);
            end
            enc_psth(1, 2:3) = enc_psth_binned(1, 2:3);
            
%             psthB = Utilities.Plotting.BinRasterSpikes(CR_psth{1, 1}, plotOptions.binSize, plotOptions.stepSize);
%             psthBS = [];
%             for row = 1:size(psthB, 1)
%                 psthBS(row, :) =  (Utilities.Smoothing.fastsmooth(psthB(row,:),plotOptions.smoothSize, 1, 0)*1e3);
%             end
%             CR_psth_binned{1, 2} = psthBS;
            
            % smooth before binning
            CR_psth_binned{1, 2} = Utilities.Plotting.BinRasterSpikes(CR_psth{1, 2}, plotOptions.binSize, plotOptions.stepSize);

            if isfield(plotOptions, 'stepSize')
                CR_psth_binned{1, 3} = CR_psth{1, 3}(1:plotOptions.stepSize:end);
            else
                CR_psth_binned{1, 3} = CR_psth{1, 3}(1:plotOptions.binSize:end);
            end
            CR_psth(1, 2:3) = CR_psth_binned(1, 2:3);
            
        end
        EncodingOrder = RecallData.EncodingOrder;
        CROrder = RecallData.CROrder;
        offsetEnc = RecallData.offsetEnc;
        offsetTones = RecallData.offsetTones;
        stimuli = RecallData.stimuli;
        
        % subbing baseline in for FR 
%         FR_psth = BTimeCourse(cellIndex, 1:3);
%         FROrder = ones(length(RecallData.stimuli), 1);
%         offsetFR = RecallData.offsetFR;        
%         f = Utilities.Plotting.PlotPerCellAllStim_Im(enc_psth, CR_psth, EncodingOrder, CROrder, offsetEnc, offsetTones, stimuli, FR_psth, FROrder, offsetFR);
%         for cl = 1:length(mcs)

        %  Pair with double mod -- current ColorMap ---------------------
        plotOptions.cols = colormap(brewermap(length(stimuli)*2, mcs{12})); % this creates a fig annoyingly
        plotOptions.cols = plotOptions.cols(2:2:end, :); 
        plotOptions.cols(1, :) = [0.4 0.4 0.4];
        plotOptions.cols(8, :) = [0.75 0 0.75];     
        %  making pair lighter and resetting grey
        plotOptions.cols = plotOptions.cols*1.45;
        plotOptions.cols(plotOptions.cols > 1) = 1; 
        plotOptions.cols(1, :) = [0.4 0.4 0.4];
        close all
        % ----------------------------------------------------------------
        
        % Accent - with replaced yellow, black and darker cols 
%         plotOptions.cols = colormap(brewermap(length(stimuli), 'Accent'));
%         plotOptions.cols(1, :) = plotOptions.cols(1, :)*0.6; % this creates a fig annoyingly        
%         plotOptions.cols(2, :) = plotOptions.cols(2, :)*0.7; % this creates a fig annoyingly        
%         plotOptions.cols(3, :) = plotOptions.cols(3, :)*0.7; % this creates a fig annoyingly        
%         plotOptions.cols(8, :) = [0.2 0.2 0.2]; % this creates a fig annoyingly        
%         plotOptions.cols(4, :) = [1 0.82 0.18]; % this creates a fig annoyingly

        % Darker set2 colors with mods - basically remakes Dark2...LOL
%         plotOptions.cols = colormap(brewermap(length(stimuli), 'Set2'));
%         plotOptions.cols = plotOptions.cols*0.85;
%         plotOptions.cols(1, :) = [0 0.6 0.6];
%         plotOptions.cols(8, :) = [0.35 0.35 0.35]; 
%         plotOptions.cols(6, :) = [1 0.82 0.18];
%         plotOptions.cols(7, :) = [0.6 0.3 0];
        
          % Linspec with mod
%         plotOptions.cols = linspecer(length(stimuli), 'qualitative'); % this creates a fig annoyingly
%         plotOptions.cols(5, :) = [0.15 0.15 0.15]; 
%         plotOptions.cols(8, :) = [0.5 0.5 0.5]; 
%         close all
        if plotOptions.showReacTrials
            [sOrd, cOrd] = sort(RecallData.CROrder); % need to sort here
            plotOptions.reacTrials = reacTrials{idx, 2}(cOrd);
        end
        f = Utilities.Plotting.PlotPerCellAllStim_Im(plotOptions, enc_psth, CR_psth, EncodingOrder, CROrder, offsetEnc, offsetTones, stimuli);
        set(findobj(f,'type','axes'),'FontName','Arial','FontSize',16,'FontWeight','Bold', 'LineWidth', 1.2);
        if plotOptions.legend
            set(findobj(f,'type','legend'),'FontName','Arial','FontSize',14,'FontWeight','Bold', 'LineWidth', 1.2);
        end
%         del(findobj(f,'type','legend')) % this suddenly doesn't work
%         sgt = sgtitle({[num2str(strctCells(cellIndex).Name) '\_' strctCells(cellIndex).brainArea '\_Comparison']}); % backslash allows you to print the underscore
%         sgt = sgtitle({[num2str(strctCells(cellIndex).Name) '\_' strctCells(cellIndex).brainArea '\_' mcs{cl}]}); % backslash allows you to print the underscore
        if plotOptions.showReacTrials
            sgt = sgtitle({[num2str(strctCells(cellIndex).Name) '\_' strctCells(cellIndex).brainArea '\_reacTrialsMarked']}); % backslash allows you to print the underscore
        else
            sgt = sgtitle({[num2str(strctCells(cellIndex).Name) '\_' strctCells(cellIndex).brainArea]}); % backslash allows you to print the underscore
        end
        sgt.FontSize = 20;
        sgt.FontWeight = 'Bold';
        
        if ~ITOnly
            pathOut = [diskPath filesep task filesep 'EncodingandCR_comparison' filesep 'AllCells'];
        else
            pathOut = [diskPath filesep task filesep 'EncodingandCR_comparison' filesep 'ITCells'];
        end
        
        if plotOptions.binSpikes
            if isfield(plotOptions, 'stepSize') && ~isequal(plotOptions.binSize, plotOptions.stepSize)
                pathOut = [pathOut '_' num2str(plotOptions.binSize) 'msbinned_' num2str(plotOptions.stepSize) 'msOverlap'];
            else
                pathOut = [pathOut '_' num2str(plotOptions.binSize) 'msbinned'];
            end
        end
%         pathOut = [diskPath filesep task filesep 'ColorTestingForAdvisors'];
%         pathOut = [diskPath filesep task filesep 'PoissonReacTesting'];
        if ~exist(pathOut, 'dir')
            mkdir(pathOut)
        end
        
        if sum(ismember(RC_names, strctCells(cellIndex).Name)) ~= 0 && isequal(strctCells(cellIndex).ChannelNumber, RC_ChanNum(find(ismember(RC_names, strctCells(cellIndex).Name))))
            pathOut = [pathOut filesep 'ReactivatedCells'];
            if ~exist(pathOut, 'dir')
                mkdir(pathOut)
            end           
        end
        if plotOptions.showReacTrials
            filename = [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name) '_ReacTrialsMarked'];
        else
            filename = [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name) '_Comparison'];
        end
%       for 173 cells took 900s at screen res and 1500s at 300dpi
%         print(f, filename, '-dsvg','-r300');
        print(f, filename, '-dpng','-r300');
%         print(f, filename, '-dpng','-r0');
        close all
%         end
    end
    
end
toc