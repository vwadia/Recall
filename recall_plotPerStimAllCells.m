
% Script to plot perStimAcrossCells
% Need to run RecallScript_main first
% currently it does this per region
% Need to change it so it loads in reacData and then plots all reac neurons
% per stim
% Include option to normalize FR
% vwadia June2023

setDiskPaths

taskCodePath = [boxPath filesep 'RecallTaskVarun'];
addpath(taskCodePath)
setTTLCodes

% set session list
task = 'Recall_Task';
sessID = Utilities.sessionListAllTasks(task, false);
sessID{2} = ['Recall_Task' filesep 'P76CS' filesep 'Recall_Session_2_20210925']; % want the rasters so using recall session
sessID = sessID';

for ss = 1%1:length(sessID)
    %% set up extra data - could save here but leave for now - June 13th 2023
    
    if ~exist('RecallData', 'var')
        load([sessID{ss} filesep 'RecallData_NoFreeRec.mat']);
    end
    tic
    % CR
    RecallData.perStimAllCellsCR = cell(length(RecallData.stimuli), 3);
    RecallData.order_perStimAllCellsCR = repelem(1:length(strctCells), RecallData.CR_numRepetitions)';
    
    rl_cr = mode(cell2mat(cellfun(@(x) length(x), RecallData.CRTimeCourse(:, 1), 'UniformOutput', false)));
    
    % per stim all cells
    for stim = 1:length(RecallData.stimuli)
        for rep = 1:length(strctCells)
            stimRaster = RecallData.CRTimeCourse{rep, 1}(find(RecallData.CROrder == stim), 1:sum(RecallData.offsetTones));
            stimPsth = RecallData.CRTimeCourse{rep, 2}(find(RecallData.CROrder == stim), 1:sum(RecallData.offsetTones));
            times = RecallData.CRTimeCourse{rep, 3};
            
            stimRaster = stimRaster(:, 1:rl_cr);
            stimPsth = stimPsth(:, 1:rl_cr);
            times = times(1:rl_cr);
            
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
    
    rl_enc = mode(cell2mat(cellfun(@(x) length(x), RecallData.EncodingTimeCourse(:, 1), 'UniformOutput', false)));
    
    
    % per stim all cells
    for stim = 1:length(RecallData.stimuli)
        for rep = 1:length(strctCells)
            stimRaster = RecallData.EncodingTimeCourse{rep, 1}(find(RecallData.EncodingOrder == stim), :);
            stimPsth = RecallData.EncodingTimeCourse{rep, 2}(find(RecallData.EncodingOrder == stim), :);
            times = RecallData.EncodingTimeCourse{rep, 3};
            
            stimRaster = stimRaster(:, 1:rl_enc);
            stimPsth = stimPsth(:, 1:rl_enc);
            times = times(1:rl_enc);
            
            %         if rep == 11 && strcmp(paths.sessPath, 'ReScreenRecall_Session_3_20210927')
            %             RecallData.perStimAllCellsEncoding{stim, 1} = [RecallData.perStimAllCellsEncoding{stim, 1}; stimRaster(:, 1:end-1)];
            %             RecallData.perStimAllCellsEncoding{stim, 2} = [RecallData.perStimAllCellsEncoding{stim, 2}; stimPsth(:, 1:end-1)];
            %         else
            %             RecallData.perStimAllCellsEncoding{stim, 1} = [RecallData.perStimAllCellsEncoding{stim, 1}; stimRaster];
            %             RecallData.perStimAllCellsEncoding{stim, 2} = [RecallData.perStimAllCellsEncoding{stim, 2}; stimPsth];
            %         end
            
            RecallData.perStimAllCellsEncoding{stim, 1} = [RecallData.perStimAllCellsEncoding{stim, 1}; stimRaster];
            RecallData.perStimAllCellsEncoding{stim, 2} = [RecallData.perStimAllCellsEncoding{stim, 2}; stimPsth];
            
            
            if rep == 1
                RecallData.perStimAllCellsEncoding{stim, 3} = [RecallData.perStimAllCellsEncoding{stim, 3}; times];
            end
        end
    end
    toc
    
    %% plotting
    
    noFR = true;
    
    if exist('allCells', 'var')
        strctCells = allCells; % resetting
    end
    
    % plot options
    perRegion = false; % to plot per region or just all cells
    reactOnly = true; % keep reactivated cells
    sepCellColor = true; % when plotting all cells together do you color all separately or color all cells from a region the same? 
    
    plotOptions.MarkerSize = 4;
    plotOptions.Fontsize = 16;
    plotOptions.LineWidth = 3;
    plotOptions.NormPsthFR = false; % normalize the FR for psth so the plots are legible
    plotOptions.visible = true;
    plotOptions.fullscreen = true;
    plotOptions.legend = true;
    
    perStimAllCellsEncoding = RecallData.perStimAllCellsEncoding;
    perStimAllCellsCR = RecallData.perStimAllCellsCR;
    
    manChoice = true; % because it is goddamnit

    if ~noFR
        perStimAllCellsFR = RecallData.perStimAllCellsFR;
    end
    
    if reactOnly
         ctype = 'ReactCells';
        if manChoice 
            
            if ss == 6
            cellstoKeep = [16; 24];
%             plotOptions.cols = [0 0.4470 0.7410; 0.4660 0.6740 0.1880];
            elseif ss == 1
                cellstoKeep = [31; 45; 106];
            end
        else
        
       
        load([diskPath filesep 'Recall_Task' filesep 'ReactivationCorrelationData0.05.mat']);
        
        % check if those cells are in here
        cellNames = cat(1, strctCells(:).Name);
        cellChans = cat(1, strctCells(:).ChannelNumber);

        reacNames = reacData.Name(:);
        reacChans = reacData.Channel(:);
        reacAreas = reacData.Area(:);
        t1 = ismember(reacNames, cellNames);
        t2 = ismember(reacChans, cellChans);

        t_ids = find(t1 & t2); % redundancy
        
        % make sure channel numbers are the same too
        for i = 1:length(t1)
          
            if t1(i) == 1
                rel_idx = find(cellNames(:, 1) == reacNames(i, 1));
                if length(rel_idx) > 1
                    if sum(reacChans(i, 1) == cellChans(rel_idx)) == 0
                        t1(i) = false;
                    end
                else
                    if ~isequal(reacChans(i, 1), cellChans(rel_idx, 1))
                        t1(i) = false;
                    end
                end
            end
        end
            
        
        t_ids = find(t1);
        assert(length(t_ids) == length(strctCells));

        relCells = reacData(t_ids, :);       
        cellstoKeep = find(reacData.reactivated(t_ids) == 1); % there are cells to keep
        
        end
        idstoKeepEnc = ismember(repelem(1:length(strctCells), RecallData.Enc_numRepetitions), cellstoKeep)';
        idstoKeepCR = ismember(repelem(1:length(strctCells), RecallData.CR_numRepetitions), cellstoKeep)';
        
        if ~noFR
            idstoKeepFR = repelem(cellstoKeep, RecallData.FR_numRepetitions);
        end
        
        allCells = strctCells;
        strctCells = strctCells(cellstoKeep);
        
        % extract relevant psths
        perStimAllCellsEncoding(:, 1:2) = cellfun(@(x) x(idstoKeepEnc, :), RecallData.perStimAllCellsEncoding(:, 1:2), 'UniformOutput', false);
            
        perStimAllCellsCR(:, 1:2) = cellfun(@(x) x(idstoKeepCR, :), RecallData.perStimAllCellsCR(:, 1:2), 'UniformOutput', false);

    else
        ctype = 'AllCells';
    end
    
    tic
   
    
    % extracting useful information from strctCells
    strctCELL = struct2cell(strctCells')';
    
    cellNames = strctCELL(:, 1); % for legend
    
    popOrder = cell2mat(strctCELL(:, 3)); % brainAreaIndices
    if ~noFR
        RecallData.RegionOrderFR = repelem(popOrder, RecallData.FR_numRepetitions);
    end
    RecallData.RegionOrderCR = repelem(popOrder, RecallData.CR_numRepetitions); % dafuq is this?
    RecallData.RegionOrderEnc = repelem(popOrder, RecallData.Enc_numRepetitions);
    
    % list of brain regions
    brainAreas = unique(popOrder, 'stable');
    regionNames = unique(strctCELL(:, 4), 'stable');
    
    % load in images
    pathFrames = [diskPath filesep sessID{ss} filesep 'stimuliUsedRecall'];
    imgs = Utilities.readInFiles(pathFrames, 'tif');
    num_imgs = length(imgs);
    ctr = 1;
    images = cell(1, length(1:num_imgs));
    for ii = 1:num_imgs
        currentfilename = [imgs(ii).folder filesep imgs(ii).name];
        currentimage = imread(currentfilename);
        images{ii} = currentimage;
        
    end
    
    stimOnCR = 5000;
    
    if noFR
        subPlotNum = 2; % 3 conditions
    else
        subPlotNum = 3; % 3 conditions
    end
    
    
    for stimIndex = l(RecallData.stimuli)
        all_colors = Utilities.distinguishable_colors(length(strctCells)); % param
        full_enc_psth = perStimAllCellsEncoding(stimIndex, 1:3);
        full_CR_psth = perStimAllCellsCR(stimIndex, 1:3);
        
        if ~noFR
            full_FR_psth = perStimAllCellsFR(stimIndex, 1:3);
        end
        
        enc_psth = cell(1, 2);
        CR_psth = cell(1, 2);
        stimuli = {};
        areaNames = {};
        if perRegion  
            for reg = 1:length(brainAreas)
                colors = all_colors(find(popOrder == brainAreas(reg)), :);
                % collect psths per region - params
               
                enc_psth = {full_enc_psth{1, 1}(find(RecallData.RegionOrderEnc == brainAreas(reg)), :), full_enc_psth{1, 2}(find(RecallData.RegionOrderEnc == brainAreas(reg)), :), full_enc_psth{1, 3}};
                CR_psth = {full_CR_psth{1, 1}(find(RecallData.RegionOrderCR == brainAreas(reg)), :), full_CR_psth{1, 2}(find(RecallData.RegionOrderCR == brainAreas(reg)), :), full_CR_psth{1, 3}};
                if ~noFR
                    FR_psth = {full_FR_psth{1, 1}(find(RecallData.RegionOrderFR == brainAreas(reg)), :), full_FR_psth{1, 2}(find(RecallData.RegionOrderFR == brainAreas(reg)), :), full_FR_psth{1, 3}};                   
                end
                
                EncodingOrder = repelem(1:sum(ismember(popOrder, brainAreas(reg))), RecallData.Enc_numRepetitions)';
                CROrder = repelem(1:sum(ismember(popOrder, brainAreas(reg))), RecallData.CR_numRepetitions)';
                offsetEnc = RecallData.offsetEnc;
                offsetTones = RecallData.offsetTones;
                areacellNames = cellNames(find(popOrder == brainAreas(reg)));
                stimuli = cellfun(@(x) num2str(x), areacellNames, 'UniformOutput', false); 
             
                % plot
                f = Utilities.Plotting.PlotPerCellAllStim_Im(plotOptions, enc_psth, CR_psth, EncodingOrder, CROrder, offsetEnc, offsetTones, stimuli);
                sgt = sgtitle({[regionNames{reg} '\_' RecallData.stimuli{stimIndex}]}); % backslash allows you to print the underscore
                sgt.FontSize = Fontsize;
                sgt.FontWeight = 'Bold';
                
                
               pathOut = [paths.basePath filesep paths.taskPath filesep paths.patientPath filesep paths.sessPath...
                    filesep paths.dataPath filesep 'PerStimPerRegionComparison'];
                if ~exist(pathOut)
                    mkdir(pathOut)
                end
                filename = [pathOut filesep [regionNames{reg} '_' RecallData.stimuli{stimIndex}]];
                
            end
        else % plot everything together
            
            % sort list of brain regions - just to have some structure
            % else it'll be ordered like strctCells - nightmare
            [blist, bid] = sort(popOrder);%, 'stable');
            regionNames = unique(strctCELL(bid, 4), 'stable');
            brainAreas = unique(popOrder); % because this is sorted can just do this
            
            for reg = 1:length(brainAreas)
                regEnc_psth = {full_enc_psth{1, 1}(find(RecallData.RegionOrderEnc == brainAreas(reg)), :), full_enc_psth{1, 2}(find(RecallData.RegionOrderEnc == brainAreas(reg)), :), full_enc_psth{1, 3}};
                regCR_psth = {full_CR_psth{1, 1}(find(RecallData.RegionOrderCR == brainAreas(reg)), :), full_CR_psth{1, 2}(find(RecallData.RegionOrderCR == brainAreas(reg)), :), full_CR_psth{1, 3}};
                if ~noFR
                    regFR_psth = {full_FR_psth{1, 1}(find(RecallData.RegionOrderFR == brainAreas(reg)), :), full_FR_psth{1, 2}(find(RecallData.RegionOrderFR == brainAreas(reg)), :), full_FR_psth{1, 3}};                   
                end
                % sort psths by region order - fix
                enc_psth{1, 1} = [enc_psth{1, 1}; regEnc_psth{1, 1}];
                enc_psth{1, 2} = [enc_psth{1, 2}; regEnc_psth{1, 2}];
                
                CR_psth{1, 1} = [CR_psth{1, 1}; regCR_psth{1, 1}];
                CR_psth{1, 2} = [CR_psth{1, 2}; regCR_psth{1, 2}];
                if reg == 1
                    enc_psth{1, 3} = regEnc_psth{1, 3};
                    CR_psth{1, 3} = regCR_psth{1, 3};
                end               
                
                % fix area cell names for legend
                areacellNames = cellNames(find(popOrder == brainAreas(reg)));
                areacellNames = cellfun(@(x) num2str(x), areacellNames, 'UniformOutput', false); 
                areaNames = [areaNames; repelem(regionNames(reg), length(areacellNames))']; 

                stimuli = [stimuli; areacellNames]; 
                
            end
            
            stimuli = cellfun(@(x, y) [y ' ' x], stimuli(:), areaNames(:),  'UniformOutput', false); 
            
            % fix order to match - if coloring all cells separately
            if sepCellColor
                EncodingOrder = repelem(1:length(strctCells), RecallData.Enc_numRepetitions)';
                CROrder = repelem(1:length(strctCells), RecallData.CR_numRepetitions)';
                pathOut = [diskPath filesep sessID{ss} filesep 'PerStim' ctype 'Comparison'];
            else
                EncodingOrder = repelem(blist, RecallData.Enc_numRepetitions)';
                CROrder = repelem(blist, RecallData.CR_numRepetitions)';
                pathOut = [diskPath filesep sessID{ss} filesep 'PerStim' ctype 'Comparison_coloredByRegion'];
            end
            offsetEnc = RecallData.offsetEnc;
            offsetTones = RecallData.offsetTones;
            % plot
            f = Utilities.Plotting.PlotPerCellAllStim_Im(plotOptions, enc_psth, CR_psth, EncodingOrder, CROrder, offsetEnc, offsetTones, stimuli);
            set(findobj(f,'type','legend'),'FontName','Arial','FontSize',14,'FontWeight','Bold', 'LineWidth', 1.2);

            sgt = sgtitle({[RecallData.stimuli{stimIndex}]}); % backslash allows you to print the underscore
            
            if reactOnly
                sgt = sgtitle({[RecallData.stimuli{stimIndex} '- Reactivated Cells']});
            else
                sgt = sgtitle({[RecallData.stimuli{stimIndex} '- All Cells']});
            end

            sgt.FontSize = plotOptions.Fontsize;
            sgt.FontWeight = 'Bold';
            
            
           
            if ~exist(pathOut)
                mkdir(pathOut)
            end
            filename = [pathOut filesep [ctype '_' RecallData.stimuli{stimIndex}]];
        end
        if manChoice
            filename = [filename '_manChoice'];
        end
        print(f, filename, '-dpng','-r0');
        close all
    end
    toc
    
end



%% manual plotting scrath code

% % 6 subplots - 1 pair for encoding, 1 for CR and 1 for FR
%                 %         f = figure; clf;
%                 f = figure('Visible', 'off'); clf;
%                 set(gcf,'Position',get(0,'Screensize')) % display fullsize on other screen
%                 sgt = sgtitle({[regionNames{reg} '\_' RecallData.stimuli{stimIndex}]}); % backslash allows you to print the underscore
%                 sgt.FontSize = Fontsize;
%                 sgt.FontWeight = 'Bold';
%                 for spN = 1:subPlotNum
%                     if spN == 1 % Encoding
%                         %                 orderToUse = find(RecallData.RegionOrderEnc == brainAreas(reg));
%                         orderToUse = repelem(1:sum(ismember(popOrder, brainAreas(reg))), RecallData.Enc_numRepetitions)';
%                         psth = enc_psth;
%                         timelimits = RecallData.offsetEnc;
%                         ttle = {'Encoding'};
%                     elseif spN == 2 % CR
%                         %                 orderToUse = find(RecallData.RegionOrderCR == brainAreas(reg));
%                         orderToUse = repelem(1:sum(ismember(popOrder, brainAreas(reg))), RecallData.CR_numRepetitions)';
%                         psth = CR_psth;
%                         timelimits = RecallData.offsetTones;
%                         ttle = {'Cued Recall'};
%                         
%                     elseif spN == 3 % FR
%                         orderToUse = repelem(1:sum(ismember(popOrder, brainAreas(reg))), RecallData.FR_numRepetitions)';
%                         psth = FR_psth;
%                         timelimits = RecallData.offsetFR;
%                         ttle = {'Free Recall'};
%                         %                 keyboard
%                     end
%                     % smoothed FR
%                     h_1(spN) = subplot(2, subPlotNum, spN);
%                     hold on
%                     imagesTOplot = unique(orderToUse);
%                     for p1 = l(imagesTOplot)
%                         Utilities.stdshade5(psth{1, 2}(find(orderToUse == imagesTOplot(p1)), :), 0.1,...
%                             colors(mod(p1, length(imagesTOplot))+1, :), psth{1, 3}, 2);
%                     end
%                     set(gca,'FontSize',Fontsize, 'FontWeight', 'bold')
%                     ylabel('Firing Rate (Hz)','FontSize',Fontsize, 'FontWeight', 'bold');
%                     title(ttle);
%                     %             xlim([psth{1, 3}(1) psth{1, 3}(end)]);
%                     xlim([round(psth{1, 3}(1)+100) round(psth{1, 3}(end)*0.95)]);
%                     yl = ylim;
%                     ylim([0 yl(2)]);
%                     plot([0 0], [0 yl(2)], '--k', 'LineWidth', LineWidth, 'HandleVisibility', 'off');
%                     
%                     if spN == 2
%                         plot([0 0], [stimOnCR stimOnCR], '--k', 'LineWidth', LineWidth, 'HandleVisibility', 'off');
%                     end
%                     if spN == subPlotNum
%                         
%                         %                 if reg == 4
%                         %                     keyboard
%                         %                 end
%                         areacellNames = cellNames(find(popOrder == brainAreas(reg)));
%                         for i = 1:length(areacellNames)
%                             areacellNames{i} = num2str(areacellNames{i});
%                         end
%                         lgnd = legend(areacellNames);
%                         if length(areacellNames) < 15
%                             lgnd.Position = [0.912673612414962,0.629427349279353,0.070312498696148,0.253723925133093];
%                         else
%                             lgnd.Position = [0.917361111949301,0.289308193532176,0.05468749916181,0.597815275517807];
%                         end
%                         %                 h_3 = subplot('Position', [0.75, 0.6, 0.25, 0.25]); % left, bottom, width, height
%                         %                 hold on
%                         %                 imshow(images{stimIndex});
%                     end
%                     % raster
%                     h_2(spN) = subplot(2, subPlotNum, spN+subPlotNum);
%                     hold on
%                     iterSize = 0;
%                     for p2 = l(imagesTOplot)
%                         iter = find(orderToUse == imagesTOplot(p2));
%                         if p2 > 1
%                             iterSize = iterSize + length(find(orderToUse == imagesTOplot(p2-1)));
%                         else
%                             iterSize = 0;
%                         end
%                         for k = 1:size(iter, 1)
%                             try
%                                 % timelimits here is in ms
%                                 plot((find(psth{1, 1}(iter(k), :)==1)+(-timelimits(1))),...
%                                     iterSize+k,'Marker','square', 'LineStyle','none','MarkerFaceColor',colors(mod(p2, length(imagesTOplot))+1, :),...
%                                     'MarkerEdgeColor','none','MarkerSize',MarkerSize)
%                                 
%                                 hold on
%                                 
%                             end
%                         end
%                     end
%                     set(gca,'FontSize',Fontsize, 'FontWeight', 'bold')
%                     ylabel('Trials (re-ordered)','FontSize',Fontsize, 'FontWeight', 'bold');
%                     %             xlim([psth{1, 3}(1) psth{1, 3}(end)]);
%                     xlim([round(psth{1, 3}(1)+100) round(psth{1, 3}(end)*0.95)]);
%                     ylim([0 size(psth{1, 1}, 1)+1]);
%                     yl2 = ylim;
%                     plot([0 0], [yl2(1) yl2(2)+1], '--k', 'LineWidth', LineWidth, 'HandleVisibility', 'off');
%                     if spN == 2
%                         plot([0 0], [stimOnCR stimOnCR], '--k', 'LineWidth', LineWidth, 'HandleVisibility', 'off');
%                     end
%                     
%                 end
