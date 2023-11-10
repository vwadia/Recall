% Script to do imagined decoding via linear model
% imaged decoding via linear model


%% Linear model

% compile response matrix for each session
% trials x cells, add a columns of ones
% regress(params, response matrix of n-1 ims)

setDiskPaths
taskPath = 'Recall_Task';

% creates feature matrices and labels
% load([diskPath filesep taskPath filesep 'FeatureMatrices_SigRampCells_avgWholePeriod_Im.mat'])
load([diskPath filesep taskPath filesep 'AllITResponses_500Stim_Im_SigRamp.mat'])
load([diskPath filesep taskPath filesep 'AllITCells_500Stim_Im_SigRamp.mat'])


load([diskPath filesep 'Recall_Task' filesep 'ReactiveITCells_alpha0.05_500Stim_Im_SigRamp.mat']);
chosenNames = cell2mat(responses(:, 3));

all_fullnames = cat(1, strctResp(:).Name);

[sessID, ~, recalledStim] = Utilities.sessionListAllTasks('Recall_Task', false, true);
sessID{2} = ['Recall_Task' filesep 'P76CS' filesep 'Recall_Session_2_20210925'];

% image features
stimDir = '500Stimuli';
load([diskPath filesep 'ObjectSpace' filesep stimDir filesep 'scores_alexnet_comparison.mat']); % will make cell array with fc6 fc6 after relu fc7 fc7 after relu
params_full = scores{1};
[coeff, score, ~, ~, explained, mu] = pca(params_full);

load([diskPath filesep 'Object_Screening' filesep 'resp_AlexNet_fc6BeforeReLu_17856ims']) 
load([diskPath filesep 'Object_Screening' filesep 'LargeImageSet_17856ims'])
proj_into_500 = (resp_wholeimage - repmat(mu, [size(resp_wholeimage, 1) 1]))*coeff; % coeff = PCs of 500 object space

reacOnly = false;
%%


if reacOnly
    chosResp = strctResp(ismember(all_fullnames, chosenNames));
else
    chosResp = strctResp;
end

fullnames = cat(1, chosResp(:).Name);

for ss = 1:length(sessID)
    
%     if ss ~= 2
        
        % create resp_mat
        load([diskPath filesep sessID{ss} filesep 'RecallData_NoFreeRec.mat'])
        strctCELL = struct2cell(strctCells');
        strctCELL = strctCELL';
        
        IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
        strctCells = strctCells(IT_Cells);
        names = cat(1, strctCells(:).Name);
        
        num_cells = sum(ismember(fullnames, names));
        
        resp_mat = [];
        
        
        
        if ss == 11
            % stupid name repetition
            nonrepids = find(ismember(fullnames, names));
            nonrepids = nonrepids(2:end);
            resp_mat = cat(2, chosResp(nonrepids).CRResp);
        elseif ss == 2
            % stupid name repetition
            nonrepids = find(ismember(fullnames, names));
            nonrepids = nonrepids(1:end-1);
            resp_mat = cat(2, chosResp(nonrepids).CRResp);
        else
            resp_mat = cat(2, chosResp(ismember(fullnames, names)).CRResp);
        end
        
        % adding column of ones
        resp_mat = [resp_mat ones(size(resp_mat, 1), 1)];
        
        all_inds = [1:size(resp_mat, 1)]';
        n_stim = length(all_inds);
        
        ndim = ceil(0.75*num_cells); % need this to be less than the number of images
        params = score(recalledStim{ss}, 1:ndim);
        
        pred_feat = [];
        for dim = 1:ndim % per dimension
            for idx = 1:length(all_inds) % per image
                
                % using leave one out
                ind_train = setdiff(all_inds, idx);
                Y_train = params(ind_train, dim);
                
                X_train = resp_mat(ind_train, :); % images x cells includes column of 1s
                
                X_test = resp_mat(idx, :);
                
                % solve for coefficients
                [b, ~, ~, ~, stats] = regress(Y_train, X_train);
                
                % predict
                pred_feat(idx, dim) = X_test*b;
                
            end
        end
        
        proj = proj_into_500(:, 1:ndim); % now the images are projected into the space built by my 500 ims
        dd = [];
        for im = 1:n_stim
            
            % target image
            target = params(im, :);
            
            % compute distances to all projected images
            dist = proj - repmat(target, [size(proj, 1), 1]);
            
            dd(im, :) = vecnorm(dist');
            % find min
            [d1(im) idx1(im)] = min(dd(im, :));
            
        end
        
        para_bpr = proj(idx1,:); % parameters of the closest match to all stim images;
        ims_bpr = im2all_woStim(:, :, idx1);  % keep images
        
        % ---------------------------------------------------------------------
        writePath = [diskPath filesep sessID{ss} filesep 'BPRecons500_woStim'];
        if reacOnly
            writePath = [writePath filesep '_reacOnly'];
        end
        if ~exist(writePath)
            mkdir(writePath)
        end
        for i = 1:n_stim
            
            im = ims_bpr(:, :, i);
            
            filename = sprintf('%04d.tif', idx1(i)); % keep the original names
            outputFileName = [writePath filesep [num2str(i) '_' filename]];
            imwrite(im, outputFileName);
            
        end
        % ---------------------------------------------------------------------
        
        ddec = [];
        
        for im = 1:n_stim
            
            % target image
            target = pred_feat(im, :);
            
            % compute distances to all projected images
            dist = proj - repmat(target, [size(proj, 1), 1]);
            
            ddec(im, :) = vecnorm(dist');
            
            % find min
            [d2(im) idx2(im)] = min(ddec(im, :));
            
        end
        
        para_recon = proj(idx2,:); % parameters of the closest match to all stim images;
        ims_decoded = im2all_woStim(:, :, idx2); % select the images;
        
        % ---------------------------------------------------------------------
        writePath = [diskPath filesep sessID{ss} filesep 'DecodedIms500_woStim'];
        if reacOnly
            writePath = [writePath filesep '_reacOnly'];
        end
        if ~exist(writePath)
            mkdir(writePath)
        end
        for i = 1:n_stim
            
            im = ims_decoded(:, :, i);
            
            filename = sprintf('%04d.tif', idx2(i)); % keep the original names
            outputFileName = [writePath filesep [num2str(i) '_' filename]];
            imwrite(im, outputFileName);
            
        end
        % ---------------------------------------------------------------------
        
        norm_dist = []; %norm(v_recon - v_original)/norm(v_bestpossrecon - v_original)
        
        % compute normalized distance
        for im = 1:n_stim
            
            v_orig = params(im, :); % parameters of original image
            
            v_recon = para_recon(im, :);
            
            v_bestPossibleRecon = para_bpr(im, :); % is this right?          
            
            norm_dist(im) = norm(v_recon - v_orig)/norm(v_bestPossibleRecon - v_orig);
            
        end
        
        clc
        sess_dist{ss} = norm_dist;
%     end
%     clc
end

%% line plot

% cols = Utilities.distinguishable_colors(length(sess_dist));
% figure; 
% hold on
% title(['Normalized Distance from ' num2str(size(proj, 1)) ' images']);
% for ss = 1:length(sess_dist)
%     
% %    plot(sess_dist{ss}, '-o', 'Color', cols(ss, :))
%    scatter(1:length(sess_dist{ss}), sess_dist{ss}', 20, cols(ss, :), 'filled')
% %    keyboard
%     
% end

%% 

all_dist = cell2mat(sess_dist);
step_size = 2;
sctns = round(prctile(all_dist, [33, 67]), 1); % 1x2 vector

f = figure; 
hold on
title({['Normalized Distance from ' num2str(size(proj, 1)) ' images'], 'Imagination'});
histogram(all_dist(all_dist <= sctns(1)), 1:step_size:sctns(1), 'EdgeColor', [1 1 1], 'FaceColor', [1 0 0], 'FaceAlpha', 1);
histogram(all_dist(all_dist > sctns(1) & all_dist <= sctns(2)), sctns(1):step_size:sctns(2),'EdgeColor', [1 1 1], 'FaceColor', [0 0 0], 'FaceAlpha', 1);
histogram(all_dist(all_dist > sctns(2)), sctns(2):step_size:ceil(max(all_dist)),'EdgeColor', [1 1 1], 'FaceColor', [0 0 1], 'FaceAlpha', 1);

% axis limits
xl = xlim;
xlim([1 xl(2)])
yl = ylim;
ylim([0 yl(2)+1])
xticks([1 20:20:140])
set(findobj(f,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1.2);

% set(gca, 'Fontsize', 12, 'FontWeight', 'bold')

filename = [diskPath filesep taskPath filesep 'forPaper' filesep 'Hist_NormDist_LargeImageSet_Im_Sigramp'];
if reacOnly
    filename = [filename '_reacOnly'];
end
% print(f, filename, '-dpng', '-r0')
