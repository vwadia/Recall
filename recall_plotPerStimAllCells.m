% extracting useful information from strctCells
strctCELL = struct2cell(strctCells');
strctCELL = strctCELL';

cellNames = strctCELL(:, 1); % for legend

popOrder = cell2mat(strctCELL(:, 3)); % brainAreaIndices
RecallData.RegionOrderFR = repelem(popOrder, RecallData.FR_numRepetitions);
RecallData.RegionOrderCR = repelem(popOrder, RecallData.CR_numRepetitions);
RecallData.RegionOrderEnc = repelem(popOrder, RecallData.Enc_numRepetitions);

brainAreas = unique(popOrder, 'stable');
regionNames = unique(strctCELL(:, 4), 'stable');

imgs = dir(fullfile(pathFrames, '*.tif'));
num_imgs = length(imgs);
ctr = 1;
images = cell(1, length(1:num_imgs));
for ii = 1:num_imgs
    currentfilename = [imgs(ii).folder filesep imgs(ii).name];
    currentimage = imread(currentfilename);
    images{ii} = currentimage;
    
end


subPlotNum = 3; % 3 conditions

% setting up viewing parameters
MarkerSize = 4;
Fontsize = 20;
LineWidth = 3;

for stimIndex = l(RecallData.stimuli)
    all_colors = Utilities.distinguishable_colors(length(strctCells)); % param

    for reg = 1:length(brainAreas)
        colors = all_colors(find(popOrder == brainAreas(reg)), :);
        % collect psths per region - params
        full_enc_psth = RecallData.perStimAllCellsEncoding(stimIndex, 1:3);
        full_CR_psth = RecallData.perStimAllCellsCR(stimIndex, 1:3);
        full_FR_psth = RecallData.perStimAllCellsFR(stimIndex, 1:3);
        
        enc_psth = {full_enc_psth{1, 1}(find(RecallData.RegionOrderEnc == brainAreas(reg)), :), full_enc_psth{1, 2}(find(RecallData.RegionOrder == brainAreas(reg)), :), full_enc_psth{1, 3}};
        CR_psth = {full_CR_psth{1, 1}(find(RecallData.RegionOrderCR == brainAreas(reg)), :), full_CR_psth{1, 2}(find(RecallData.RegionOrder == brainAreas(reg)), :), full_CR_psth{1, 3}};
        FR_psth = {full_FR_psth{1, 1}(find(RecallData.RegionOrderFR == brainAreas(reg)), :), full_FR_psth{1, 2}(find(RecallData.RegionOrder == brainAreas(reg)), :), full_FR_psth{1, 3}};
        
        % 6 subplots - 1 pair for encoding, 1 for CR and 1 for FR
        f = figure; clf;
        set(gcf,'Position',get(0,'Screensize')) % display fullsize on other screen
        sgt = sgtitle({[regionNames{reg} '\_' RecallData.stimuli{stimIndex}]}); % backslash allows you to print the underscore
        sgt.FontSize = Fontsize;
        sgt.FontWeight = 'Bold';
        for spN = 1:subPlotNum
            if spN == 1 % Encoding
                orderToUse = find(RecallData.RegionOrderEnc == brainAreas(reg));
                psth = enc_psth;
                timelimits = RecallData.offsetEnc;
                ttle = {'Encoding'};
            elseif spN == 2 % CR
                orderToUse = find(RecallData.RegionOrderCR == brainAreas(reg));
                psth = CR_psth;
                timelimits = RecallData.offsetTones;
                ttle = {'Cued Recall'};
                
            elseif spN == 3 % FR
                orderToUse = find(RecallData.RegionOrderFR == brainAreas(reg));
                psth = FR_psth;
                timelimits = RecallData.offsetFR;
                ttle = {'Free Recall'};
%                 keyboard
            end
            % smoothed FR
            h_1(spN) = subplot(2, subPlotNum, spN);
            hold on
            imagesTOplot = unique(orderToUse);
            for p1 = l(imagesTOplot)
                Utilities.stdshade5(psth{1, 2}(find(orderToUse == imagesTOplot(p1)), :), 0.1,...
                    colors(mod(p1, length(imagesTOplot))+1, :), psth{1, 3}, 2);
            end
            set(gca,'FontSize',Fontsize, 'FontWeight', 'bold')
            ylabel('Firing Rate (Hz)','FontSize',Fontsize, 'FontWeight', 'bold');
            title(ttle);
%             xlim([psth{1, 3}(1) psth{1, 3}(end)]);
            xlim([round(psth{1, 3}(1)+100) round(psth{1, 3}(end)*0.95)]);
            yl = ylim;
            ylim([0 yl(2)]);
            plot([0 0], [0 yl(2)], '--k', 'LineWidth', LineWidth, 'HandleVisibility', 'off');
            if spN == subPlotNum
%                 if reg == 4
%                     keyboard
%                 end
                areacellNames = cellNames(find(popOrder == brainAreas(reg)));
                for i = 1:length(areacellNames)
                    areacellNames{i} = num2str(areacellNames{i});
                end
                lgnd = legend(areacellNames);
                if length(areacellNames) < 15
                    lgnd.Position = [0.912673612414962,0.629427349279353,0.070312498696148,0.253723925133093];
                else
                    lgnd.Position = [0.917361111949301,0.289308193532176,0.05468749916181,0.597815275517807];
                end
%                 h_3 = subplot('Position', [0.75, 0.6, 0.25, 0.25]); % left, bottom, width, height
%                 hold on
%                 imshow(images{stimIndex});
            end
            % raster
            h_2(spN) = subplot(2, subPlotNum, spN+subPlotNum);
            hold on
            iterSize = 0;
            for p2 = l(imagesTOplot)
                iter = find(orderToUse == imagesTOplot(p2));
                if p2 > 1
                    iterSize = iterSize + length(find(orderToUse == imagesTOplot(p2-1)));
                else
                    iterSize = 0;
                end
                for k = 1:size(iter, 1)
                    try
                        % timelimits here is in ms
                        plot((find(psth{1, 1}(iter(k), :)==1)+(-timelimits(1))),...
                            iterSize+k,'Marker','square', 'LineStyle','none','MarkerFaceColor',colors(mod(p2, length(imagesTOplot))+1, :),...
                            'MarkerEdgeColor','none','MarkerSize',MarkerSize)
                        
                        hold on
                        
                    end
                end
            end
            set(gca,'FontSize',Fontsize, 'FontWeight', 'bold')
            ylabel('Trials (re-ordered)','FontSize',Fontsize, 'FontWeight', 'bold');
%             xlim([psth{1, 3}(1) psth{1, 3}(end)]);
            xlim([round(psth{1, 3}(1)+100) round(psth{1, 3}(end)*0.95)]);
            ylim([0 size(psth{1, 1}, 1)+1]); 
            yl2 = ylim;
            plot([0 0], [yl2(1) yl2(2)+1], '--k', 'LineWidth', LineWidth, 'HandleVisibility', 'off');
            
            
        end
%         linkaxes(h_1, 'y');
        pathOut = [paths.basePath filesep paths.taskPath filesep paths.patientPath filesep paths.sessPath...
            filesep paths.dataPath filesep 'PerStimPerRegionComparison'];
        if ~exist(pathOut)
            mkdir(pathOut)
        end
        filename = [pathOut filesep [regionNames{reg} '_' RecallData.stimuli{stimIndex}]];
        print(f, filename, '-dpng','-r0');
        close all
    end
end