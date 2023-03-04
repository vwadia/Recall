% Correlate vector neurons x timebin after viewin vs
% (vector neurons x timebin after CR/FR event) <-- slide this across whole CR/FR period
% Can try with both normalized and non-normalized FRs


% This script needs: 

%% Load in the data

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

% IT cells only
load([diskPath filesep 'Recall_Task' filesep 'ITCells_500stim_Im']); % creates psths, responses, strctCells for IT cells across sessions
load([diskPath filesep 'Recall_Task' filesep 'ITResponses_500stim_Im']); % creates strctResp with screening, encoding, and CRResp for imagined images

% IT + MTL cells

% IT + MFC + MTL

% ALL CELLS
% load([diskPath filesep 'Recall_Task' filesep 'AllCells_500stim_Im']);
% load([diskPath filesep 'Recall_Task' filesep 'AllResponses_500stim_Im']);


% params
load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'params_Alexnet_fc6_500Stimuli.mat']); % will create params = 500x50

%% isolate stimuli in the 3 conditions

% strctCELLS = struct2cell(strctCells);
% strctCELLS = reshape(strctCELLS, [numel(fieldnames(strctCells)) size(strctCells, 2)]);
enc_m = {};
im_m = {};
scrn_m = {};
c_s2 = 1;
c_s3 = 1;
for cellIndex = l(strctCells)
    if strctCells(cellIndex).SessionID == 1
        enc_m{1}(cellIndex, :) = strctResp(cellIndex).EncResp;
        im_m{1}(cellIndex, :) = strctResp(cellIndex).CRResp;
        scrn_m{1}(cellIndex, :) = strctResp(cellIndex).ScrnResp;
    elseif strctCells(cellIndex).SessionID == 2
        enc_m{2}(c_s2, :) = strctResp(cellIndex).EncResp;
        im_m{2}(c_s2, :) = strctResp(cellIndex).CRResp;
        scrn_m{2}(c_s2, :) = strctResp(cellIndex).ScrnResp;
        c_s2 = c_s2+1;
    elseif strctCells(cellIndex).SessionID == 3
        enc_m{3}(c_s3, :) = strctResp(cellIndex).EncResp;
        im_m{3}(c_s3, :) = strctResp(cellIndex).CRResp;
        scrn_m{3}(c_s3, :) = strctResp(cellIndex).ScrnResp;
        c_s3 = c_s3+1;
    end
end

%% Plotting distance in neural space

% Matrix for viewed stim (cells x images) and a matrix for imagined stim
% for each column of viewed matrix (individual stim in neural space)
% compute distance between it and all columns of imagined matrix
% plot the distance matrix per session
% reinstatement score: per viewed image if the imagined image is closest =
% 1, otherwise 0. Sum for all images

dist_mat = {};

% normalize
for sess = 1:3
    % note here each row is a cell
%     enc_m{sess} = zscore(enc_m{sess}, [], 2);
%     im_m{sess} = zscore(im_m{sess}, [], 2);
%     scrn_m{sess} = zscore(scrn_m{sess}, [], 2);
    
    e_m = zscore(enc_m{sess}, [], 2);
    i_m = zscore(im_m{sess}, [], 2);
    s_m = zscore(scrn_m{sess}, [], 2);
    
    for s_1 = 1:size(i_m, 2)
        for s_2 = 1:size(i_m, 2)
%             dist_mat{sess}(s_1, s_2) = norm(e_m(:, s_1) - i_m(:, s_2));
            dist_mat{sess}(s_1, s_2) = norm(s_m(:, s_1) - i_m(:, s_2));
%             dist_mat{sess}(s_1, s_2) = norm(s_m(:, s_1) - e_m(:, s_2)); % sanity check
            
        end
        dist_mat{sess}(s_1, :) = (dist_mat{sess}(s_1, :))./max(dist_mat{sess}(s_1, :));
    end
    
    % shouldn't these matrices be symmetrical? - nope
    % 2, 1 is dist from 2nd viewed to 1st imagined
    % 1, 2 is dist from 1st viewed to 2nd imagined
    f = figure;
    set(gcf,'Position',get(0,'Screensize')) % display fullsize on other screen
    imagesc(dist_mat{sess}); colormap(flipud(parula)); 
    cb = colorbar;
    cb.Ticks = [min(min(dist_mat{sess})) max(max(dist_mat{sess}))];
    cb.Label.String = 'Normalized Distance';
    xlabel('Imagined stim');
    title('Distance between Screening And Imagined images', ['Session\_' num2str(sess)]);
    ylabel('Viewed stim (scrn)');
%     ylabel('Encoded stim');
    filename = [taskPath filesep 'DistMat_ScrnIm_P76CS_Sess_' num2str(sess)];
    print(f, filename, '-dpng', '-r0');
    close all
end


%% Nearest neighbour decoding

% Is the nearest neighbour of a viewed image it's imagined counterpart?
% add an increasing number of distractors and find min distance like Steven did in cell paper
% make this a function?

n_repeats = 1000;
n_dist = 1;
dec_acc = {};
toPlot = 0;
for sess = 1:3
    
    e_m = zscore(enc_m{sess}, [], 2);
    i_m = zscore(im_m{sess}, [], 2);
    s_m = zscore(scrn_m{sess}, [], 2);

    if toPlot
        if sess == 1 || sess == 2
            stimuli = ["1", "2", "3", "4", "5", "6",  "Chance"];
        elseif sess == 3
            stimuli = ["1", "2", "3", "4", "5", "6", "7", "8", "Chance"];
        end
        stimuli = cellstr(stimuli);
        f = figure;
        set(gcf,'Position',get(0,'Screensize')) % display fullsize on other screen
        hold on
    end
    
    for stim = 1:size(i_m, 2)
        chance = [];
        for n_dist = 1:size(i_m, 2)-1  
            chance(end+1) = 1/(n_dist+1);
            dec_acc_im = zeros(1, n_repeats);
            for rep = 1:n_repeats
                % comparing each visual image to the imagined ones
                sample_set = setdiff(1:size(i_m, 2), stim);
                sample_ids = [randsample(sample_set, n_dist, false) stim];
                
                sample_ims = i_m(:, sample_ids);
                
%                 if sess == 1 || sess == 3
                    target_im = s_m(:, stim);
%                 elseif sess == 2
%                     target_im = e_m(:, stim);
%                 end
                
                dist = [];
                for run = 1:length(sample_ids)
                    dist(run) = norm(sample_ims(:, run) - target_im);
                end
                
                % find minimum
                [min_dist, min_idx] = min(dist);
                if sample_ids(min_idx) == stim
                    dec_acc_im(rep) = 1;
                end
            end
            dec_acc{sess}(stim, n_dist) = sum(dec_acc_im(:))/n_repeats;           
        end
        if toPlot
            plot(dec_acc{sess}(stim, :), 'LineWidth', 1.5);
        end
    end
%     dec_acc{sess}(end+1, :) = chance;
    if toPlot
        plot(chance, '-k', 'LineWidth', 2);
        xlim([1 size(i_m, 2)]);
        xlabel('Number of Distracor images');
        ylabel('Decoding accuracy');
        title('NN Decoding comparing Screening And Encoding', ['Session\_' num2str(sess)]); legend(stimuli);
        filename = [taskPath filesep 'NNDec_ScrnEnc_P76CS_Sess' num2str(sess)];
        print(f, filename, '-dpng', '-r0');
        close all
    end
  
end
%% plot heatmaps for decoding

for sess = 1:3
    f = figure; 
    set(gcf,'Position',get(0,'Screensize')) % display fullsize on other screen
    imagesc(dec_acc{sess}); 
    cb = colorbar;
    cb.Label.String = 'NN Decoding Accuracy';
    cb.Ticks = [0 max(max(dec_acc{sess}))];
    
    xlabel('Number of Distracor images');
    ylabel('Stimulus Number');
    if sess == 1 || sess == 3
        title('NN Decoding comparing Screening And Imagination', ['Session\_' num2str(sess)]); 
        filename = [taskPath filesep 'NNDec_Heatmap_ScrnIm_P76CS_Sess' num2str(sess)];
    elseif sess == 2
        title('NN Decoding comparing Encoding And Imagination', ['Session\_' num2str(sess)]); 
        filename = [taskPath filesep 'NNDec_Heatmap_EncIm_P76CS_Sess' num2str(sess)];
    end
    print(f, filename, '-dpng', '-r0');
    close all

end










