% Script to plot SFC results either total summary
% or per session
% Didn't want too much clutter in the compute scripts
% note these is written for the new pipeline
% vwadia Feb2023

%% summary figure

setDiskPaths
task = 'Recall_Task';
% scale = 'ppc_linear';
scale = 'ppc_log';

HPC = false; balanced = true;
addPVals = true;
if HPC && balanced
        datPath = [diskPath filesep task filesep scale filesep 'Data' filesep 'HPCOutput' filesep 'combined'];
%         fn = 'ppc_LFFACell_LHLFP_sigRamp_ScreeningImagination_Cell_combined';
%         fn = 'ppc_LHCell_LFFALFP_sigRamp_ScreeningImagination_Cell_combined';
        fn = 'ppc_RFFACell_RHLFP_sigRamp_ScreeningImagination_Cell_combined';
%         fn = 'ppc_RHCell_RFFALFP_sigRamp_ScreeningImagination_Cell_combined';

%         fn = 'ppc_RFFACell_ROFLFP_sigRamp_ScreeningImagination_Cell_combined';
%         fn = 'ppc_RFFACell_RSMALFP_sigRamp_ScreeningImagination_Cell_combined';
%         fn = 'ppc_RFFACell_RALFP_sigRamp_ScreeningImagination_Cell_combined';

%         fn = 'ppc_LFFACell_LOFLFP_sigRamp_ScreeningImagination_Cell_combined';
%         fn = 'ppc_LFFACell_LSMALFP_sigRamp_ScreeningImagination_Cell_combined';
%         fn = 'ppc_LFFACell_LALFP_sigRamp_ScreeningImagination_Cell_combined';

        load([datPath filesep fn]);
        
        for nr = 1:size(ppc, 1)
            if ~isempty(ppc{nr, 1})
                if ndims(ppc{nr, 1}) > ndims(ppc{nr, 2})
                    toAverage = 1;
                else
                    toAverage = 2;
                end
            
                ppc{nr, toAverage} = mean(ppc{nr, toAverage}, 4);
            end
        end
else
    if strcmp(scale, 'ppc_linear')
        frq = linspace(2, 100, 98);
        load([diskPath filesep task filesep scale filesep 'ppc_RITCellRHippLFP_allCells.mat']); noSess2 = 1; Side = 'RFFA';
        
    elseif strcmp(scale, 'ppc_log')
%         datPath = [diskPath filesep task filesep scale];
        datPath = [diskPath filesep task filesep scale filesep 'SpikesfromstimOn'];

        
        %     load([datPath filesep 'ppc_RFFACell_RHLFP_sigRamp_EncodingImagination_cellChans']); Side = 'RFFA';
        load([datPath filesep 'ppc_RFFACell_RHLFP_sigRamp_ScreeningImagination_cellChans']); Side = 'RFFA';
%         load([datPath filesep 'ppc_RHCell_RFFALFP_sigRamp_ScreeningImagination_cellChans']); Side = 'RFFA';
        
    end
end

cols = Utilities.distinguishable_colors(2);
avPPC_Sc = [];
avPPC_Im = [];
ppc_Sc = ppc(:, 1);
ppc_Im = ppc(:, 2);
for s = 1:length(ppc_Sc)
    if ~isempty(ppc_Sc{s}) && ~isempty(ppc_Im{s})
        if size(ppc_Sc{s}, 1) == 1 && size(ppc_Sc{s}, 2) == 1
            avPPC_Sc = [avPPC_Sc squeeze(ppc_Sc{s})];
        elseif size(ppc_Im{s}, 1) == 1 && size(ppc_Im{s}, 2) == 1
            avPPC_Im = [avPPC_Im squeeze(ppc_Im{s})];
        elseif size(ppc_Sc{s}, 1) == 1 || size(ppc_Sc{s}, 2) == 1
            avPPC_Sc = [avPPC_Sc squeeze(nanmean(ppc_Sc{s}))];
        elseif size(ppc_Im{s}, 1) == 1 || size(ppc_Im{s}, 2) == 1
            avPPC_Im = [avPPC_Im squeeze(nanmean(ppc_Im{s}))];
        else
            avPPC_Sc = [avPPC_Sc squeeze(nanmean(nanmean(ppc_Sc{s})))];
            avPPC_Im = [avPPC_Im squeeze(nanmean(nanmean(ppc_Im{s})))];
        end                
    end
end

% old way of plotting std error of nanmean - just use stdshade instead ------
% avPPC_Sc = nanmean(avPPC_Sc, 2);
% avPPC_Im = nanmean(avPPC_Im, 2);
% stdErr_Sc = std(avPPC_Sc)/sqrt(length(avPPC_Sc));
% stdErr_Im = std(avPPC_Im)/sqrt(length(avPPC_Im));
% ------------------------------------------------------------------------


if addPVals
    
    ppc_C1_all = [];
    ppc_C2_all = [];
    for session = 1:length(ppc)
        
        if ~isempty(ppc{session})
            
            ppc_C1 = ppc{session, 1};
            n_neurons_C1 = size(ppc_C1, 1);
            n_channels_C1 = size(ppc_C1, 2);
            n_freq_C1 = size(ppc_C1, 3);
            
            ppc_C1 = reshape(ppc_C1, [n_neurons_C1*n_channels_C1 n_freq_C1]);
            
            ppc_C2 = ppc{session, 2};
            n_neurons_C2 = size(ppc_C2, 1);
            n_channels_C2 = size(ppc_C2, 2);
            n_freq_C2 = size(ppc_C2, 3);
            
            ppc_C2 = reshape(ppc_C2, [n_neurons_C2*n_channels_C2 n_freq_C2]);
            
            frqs = frq;
            ppc_C1_all = [ppc_C1_all; ppc_C1];
            ppc_C2_all = [ppc_C2_all; ppc_C2];
        end
        
    end
    
    for i = 1:n_freq_C1
        
        [~, pvals(i)] = ttest(ppc_C1_all(:, i), ppc_C2_all(:, i));
        
    end
    
    sumFig = figure;
    h1 = subplot(2, 1, 1);
    
    hold on
    Utilities.stdshade5(avPPC_Sc', 0.1, cols(1, :), frq');
    Utilities.stdshade5(avPPC_Im', 0.1, cols(2, :), frq');
%     xlim([params.low_freq 40])
    set(gca, 'FontSize', 9, 'FontWeight', 'bold')

    h2 = subplot(2, 1, 2);
%     hold on
    semilogy(frq, pvals, 'Color', [0 0.4 0], 'LineWidth', 2);
    yline(0.05/97, 'LineWidth', 1.5)
    ylim([0 0.9])
%     xlim([params.low_freq 40])

    linkaxes([h1, h2], 'x')
    sgt = sgtitle([params.cellArea ' neuron to ' params.lfpArea ' micro-lfp']);
    sgt.FontWeight = 'bold';

    set(gca, 'FontSize', 9, 'FontWeight', 'bold')

    filename = [datPath filesep [params.cellArea 'Cell_' params.lfpArea '_micro-lfp_wPVals']];
    print(sumFig, filename, '-dpng', '-r0')

else
sumFig = figure;
% plot(f_Sc{1}, avPPC_Sc, 'LineWidth', 2)
% plot(f_Im{1}, avPPC_Im, 'LineWidth', 2)
% if strcmp(scale, 'ppc_linear')
    hold on
    Utilities.stdshade5(avPPC_Sc', 0.1, cols(1, :), frq');
    Utilities.stdshade5(avPPC_Im', 0.1, cols(2, :), frq');
% elseif strcmp(scale, 'ppc_log')
%     % plotting this with shaded error is a pain
%     s1 = semilogx(frq, nanmean(avPPC_Sc, 2), frq, nanmean(avPPC_Im, 2), 'LineWidth', 2);
% end
    title([params.cellArea ' neuron to ' params.lfpArea ' micro-lfp'])
%     ylim([0 10e-4])
% xlim([params.low_freq 40])
filename = [datPath filesep [params.cellArea 'Cell_' params.lfpArea '_micro-lfp']];
% filename = [diskPath filesep task filesep 'PPC_' Side 'Cell_' Side(1) 'HLFP_AllSessions_AllResp'];
legend('BottomUp', 'TopDown');
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
% print(sumFig, filename, '-dpng', '-r0')

end
%% per session
setDiskPaths
task = 'Recall_Task';

% scale = 'ppc_linear';
scale = 'ppc_log';

% if strcmp(scale, 'ppc_linear')
%     frq = linspace(2, 100, 98);
%     load([diskPath filesep task filesep scale filesep 'ppc_RITCellRHippLFP_allCells.mat']); noSess2 = 1; Side = 'RFFA';
% 
% elseif strcmp(scale, 'ppc_log')
% %     load([diskPath filesep task filesep scale filesep 'ppc_RFFACell_RHLFP_sigRamp_EncodingImagination_cellChans']); Side = 'RFFA';
%     load([diskPath filesep task filesep scale filesep 'ppc_RFFACell_RHLFP_sigRamp_ScreeningImagination_cellChans']); Side = 'RFFA';
%     
% end


cols = Utilities.distinguishable_colors(2);
avPPC_Sc = [];
avPPC_Im = [];
for session = 1:length(ppc_Sc)
    if session == 7
        keyboard
    end
    if ~isempty(ppc_Sc{session}) && ~isempty(ppc_Im{session})
        if size(ppc_Sc{session}, 1) == 1 || size(ppc_Sc{session}, 2) == 1
            avPPC_Sc = squeeze(ppc_Sc{session});        
        elseif size(ppc_Sc{session}, 1) == 1 && size(ppc_Sc{session}, 2) == 1
            avPPC_Sc = reshape(ppc_Sc{session}, [length(ppc_Sc{session}) 1]);  
        else
            avPPC_Sc = squeeze(nanmean(ppc_Sc{session}));
        end
        
        if size(ppc_Im{session}, 1) == 1 || size(ppc_Im{session}, 2) == 1
            avPPC_Im = squeeze(ppc_Im{session});
        elseif size(ppc_Im{session}, 1) == 1 && size(ppc_Im{session}, 2) == 1
            avPPC_Im = reshape(ppc_Im{session}, [length(ppc_Im{session}) 1]);
        else
            avPPC_Im = squeeze(nanmean(ppc_Im{session}));
        end
        f = figure; 
%         if strcmp(scale, 'ppc_linear')
            hold on

            Utilities.stdshade5(avPPC_Sc, 0.1, cols(1, :), frq');
            Utilities.stdshade5(avPPC_Im, 0.1, cols(2, :), frq');
%         elseif strcmp(scale, 'ppc_log')
%            c1_ppc = nanmean(avPPC_Sc);
%            stderr_c1 = nanmean(c1_ppc)/sqrt(length(c1_ppc));
%            c2_ppc = nanmean(avPPC_Im);
%            stderr_c2 = nanmean(c2_ppc)/sqrt(length(c2_ppc));
%            semilogx(frq, c1_ppc, frq, c2_ppc, 'LineWidth', 2)
% %            patch('Faces', [frq fliplr(frq)], [c1_ppc + stderr_c1 fliplr(c1_ppc - stderr_c1)], 'FaceColor', [0 0 1],'FaceAlpha', 0.1)
%         end
        
%         xticks(frq);
%         set(gca, 'xticklabel', {[]}); % clear the old labels
%         xtik = get(gca,'xtick');
%         xtiklabs = {frq};
        yl = ylim;
        ylim([0 yl(2)])
        if exist('cellChans', 'var')
            filename = [diskPath filesep task filesep scale filesep 'ppc_' Side 'Cell_' Side(1) 'HLFP_Sess' num2str(session) '_allCells_cellChans'];
%             filename = [diskPath filesep task filesep scale filesep 'ppc_' Side 'Cell_' Side(1) 'HLFP_Sess' num2str(session) '_sigRamp_cellChans'];            
        else
            filename = [diskPath filesep task filesep scale filesep 'ppc_' Side 'Cell_' Side(1) 'HLFP_Sess' num2str(session) '_allCells'];
%             filename = [diskPath filesep task filesep scale filesep 'ppc_' Side 'Cell_' Side(1) 'HLFP_Sess' num2str(session) '_sigRamp'];
        end
        legend('BottomUp', 'TopDown');
        set(gca, 'FontSize', 14, 'FontWeight', 'bold')
%         print(f, filename, '-dpng', '-r0')
%         close all
    end
    
end