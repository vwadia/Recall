
%% MAKE THIS A PLOTTING FUNCTION

% Plot 1
% per cell encoding, CR and FR on the same page for all stimuli
noFR = true;

if noFR
    subPlotNum = 2; % only enc and CR
else
    subPlotNum = 3;
end

% setting up viewing parameters
MarkerSize = 4;
Fontsize = 20;
LineWidth = 3;

for cellIndex = l(strctCells)
    
    colors = Utilities.linspecer(length(RecallData.stimuli));
    
    for stim = 1:length(RecallData.stimuli)
        % collect psths
        full_enc_psth = RecallData.EncodingTimeCourse(cellIndex, 1:3);
        full_CR_psth = RecallData.CRTimeCourse(cellIndex, 1:3);
        
        if ~noFR
            full_FR_psth = RecallData.perCellAllStimFR(cellIndex, 1:3);
        end
        
        enc_psth = {full_enc_psth{1, 1}(find(RecallData.EncodingOrder == stim), :), full_enc_psth{1, 2}(find(RecallData.EncodingOrder == stim), :), full_enc_psth{1, 3}};
        CR_psth = {full_CR_psth{1, 1}(find(RecallData.CROrder == stim), :), full_CR_psth{1, 2}(find(RecallData.CROrder == stim), :), full_CR_psth{1, 3}};
        if ~noFR
            FR_psth = {full_FR_psth{1, 1}(find(RecallData.order_perCellAllStimFR == stim), :), full_FR_psth{1, 2}(find(RecallData.order_perCellAllStimFR == stim), :), full_FR_psth{1, 3}};
        end
        
        % 6 subplots - 1 pair for encoding, 1 for CR and 1 for FR
        f = figure; clf;
        set(gcf,'Position',get(0,'Screensize')) % display fullsize on other screen
        sgt = sgtitle({[num2str(strctCells(cellIndex).Name) '\_' strctCells(cellIndex).brainArea '\_' RecallData.stimuli{stim}]}); % backslash allows you to print the underscore
        sgt.FontSize = Fontsize;
        sgt.FontWeight = 'Bold';
        
        for spN = 1:subPlotNum
            if spN == 1 % Encoding
                orderToUse = repelem(stim, size(enc_psth{1, 1}, 1))';
%                 psth = {full_enc_psth{1, 1}(find(RecallData.EncodingOrder == stim), :), full_enc_psth{1, 2}(find(RecallData.EncodingOrder == stim), :), full_enc_psth{1, 3}};
                psth = enc_psth;
                timelimits = RecallData.offsetEnc;
                ttle = {'Encoding'};
            elseif spN == 2 % CR
                orderToUse = repelem(stim, size(CR_psth{1, 1}, 1))';
%                 psth = {full_CR_psth{1, 1}(find(RecallData.CROrder == stim), :), full_CR_psth{1, 2}(find(RecallData.CROrder == stim), :), full_CR_psth{1, 3}};
                psth = CR_psth;
                timelimits = RecallData.offsetTones;
                ttle = {'Cued Recall'};
            elseif spN == 3 % FR
                orderToUse = repelem(stim, size(FR_psth{1, 1}, 1))';
%                 psth = {full_FR_psth{1, 1}(find(RecallData.order_perCellAllStimFR == stim), :), full_FR_psth{1, 2}(find(RecallData.order_perCellAllStimFR == stim), :), full_FR_psth{1, 3}};
                psth = FR_psth;
                timelimits = RecallData.offsetFR;
                ttle = {'Free Recall'};
            end
            % smoothed FR
            h_1(spN) = subplot(2, subPlotNum, spN);
            hold on
            imagesTOplot = unique(orderToUse);
            
            if noFR
                gyl = Utilities.Plotting.findingGlobalYLim(CR_psth{1, 2}(:, 1:6000), 2, repelem(2, 6), 'AIC');
            else
                gyl = Utilities.Plotting.findingGlobalYLim([CR_psth{1, 2}(:, 1:6000); FR_psth{1, 2}(:, 1:6000)], [2; 3], repelem([2; 3], 6), 'AIC');
            end
            
            gyl2 = Utilities.Plotting.findingGlobalYLim(enc_psth{1, 2}, 1, repelem(1, 6)', 'AIC');
            % set globaly yl
            if gyl2 >= gyl
                globalyl = gyl2;
            else
                globalyl = gyl;
            end
            % set globaly yl
            if globalyl == 0
                globalyl = 1;
            end
            
            for p1 = l(imagesTOplot)
                Utilities.stdshade5(psth{1, 2}(find(orderToUse == imagesTOplot(p1)), :), 0.1,...
                    colors(stim, :), psth{1, 3}, 2);
            end
            set(gca,'FontSize',Fontsize, 'FontWeight', 'bold')
            ylabel('Firing Rate (Hz)','FontSize',Fontsize, 'FontWeight', 'bold');
            title(ttle);
            xlim([round(psth{1, 3}(1)+100) round(psth{1, 3}(end)*0.95)]);
            yl = ylim;
            ylim([0 globalyl]);
            plot([0 0], [0 globalyl], '--k', 'LineWidth', LineWidth, 'HandleVisibility', 'off');
            if spN == subPlotNum
                lgnd = legend(RecallData.stimuli{stim});
                lgnd.Position = [0.912673612414962,0.629427349279353,0.070312498696148,0.253723925133093];
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
                            iterSize+k,'Marker','|', 'LineStyle','none','MarkerFaceColor',colors(stim, :),...
                            'MarkerEdgeColor',colors(stim, :),'MarkerSize',24,  'linewidth', 1.5)
                        
                        hold on
                        
                    end
                end
            end
            set(gca,'FontSize',Fontsize, 'FontWeight', 'bold')
            ylabel('Trials (re-ordered)','FontSize',Fontsize, 'FontWeight', 'bold');
            xlim([round(psth{1, 3}(1)+100) round(psth{1, 3}(end)*0.95)]);
            ylim([0 size(psth{1, 1}, 1)+1]);
            plot([0 0], [0 size(psth{1, 1}, 1)+1], '--k', 'LineWidth', LineWidth, 'HandleVisibility', 'off');
            
            linkaxes(h_1, 'y');
%             linkaxes(h_2, 'y');
            linkaxes([h_1(spN) h_2(spN)], 'x');
        end
        
        pathOut = [paths.basePath filesep paths.taskPath filesep paths.patientPath filesep paths.sessPath...
            filesep paths.dataPath filesep 'PerCellAllIndividualStim'];
        if ~exist(pathOut)
            mkdir(pathOut)
        end
        filename = [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name) '_' RecallData.stimuli{stim}];
        print(f, filename, '-dpng','-r0');
        close all
    end
end