
% Script for decoding imagined stimuli

% Note that this is essentialy exactly 

% Take screening axis for cell
% Project imagined responses onto it 
% decode features

setDiskPaths

taskPath = 'Object_Screening';

% load in params
layermat = 'fc6';
stimDir = '500Stimuli';
ndim = 5;
load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_Alexnet_' layermat '_' stimDir '.mat']]);

% all screening neurons (rescreen or the morning session for P76)
load([diskPath filesep 'Recall_Task' filesep 'AllITCells_500Stim_Im_SigRamp.mat'])
load([diskPath filesep 'Recall_Task' filesep 'AllITResponses_500Stim_Im_SigRamp.mat'])
fullnames = cat(1, strctResp(:).Name);

load([diskPath filesep 'Recall_Task' filesep 'FeatureMatrices_SigRampCells_avgWholePeriod_Im'])


sessID = {['Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917'],...
    ['Recall_Task' filesep 'P76CS' filesep 'Recall_Session_2_20210925'],... % want the rasters so using recall session
    ['Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_3_20210927'],...
    ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_1_20220330'],...
    ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_2_20220403'],...
    ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_3_20220405'],...
    ['Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_1_20220728'],...
    ['Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_2_20220731']};



recalledStim = {[12,19,25,123,270,487],...
    [54,129,130,186,270,449],...
    [18,44,45,81,135,181,230,344],...
    [9,157,167,200,201,291,422,498],...
    [9,12,117,292,360,368,421,492],...
    [77,112,160,232,278,345,387,440],...
    [17,61,76,114,157,161,177,480],...
    [55,88,148,251,256,274,285,365]};


%%

sess_acc = [];

for ss = 1:length(featureMatrices)
    
    tic
    % create resp_mat
    load([diskPath filesep sessID{ss} filesep 'RecallData_NoFreeRec.mat'])
    strctCELL = struct2cell(strctCells');
    strctCELL = strctCELL';
    
    IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
    strctCells = strctCells(IT_Cells);
    names = cat(1, strctCells(:).Name);
    
    valid_cells = ismember(fullnames, names);
    fr_mat = cat(2, strctResp(valid_cells).CRResp);
    
    % make sure the params are of the right dimensions
    params = params(:, 1:ndim);
    
    % zscoring features
    params = zscore(params);
    pred_feat = [];
    obs_feat = [];
    % get predicted features for each image per neuron
    for dimIndex = 1:ndim
        imageIDs = 1:length(recalledStim{ss});
        % observed responses - note zscore
        feat = params(recalledStim{ss}, dimIndex);
        
        % predicted features
        % currently I am not normalizing each img fr vector per cell...is
        % this an issue?
        % method = 'linear_regression';
        method = 'sta';
        [pred_feat(:, dimIndex), obs_feat(:, dimIndex)] = Utilities.computePredictedNNFeatures(feat, fr_mat, imageIDs, method);
        
    end
    
    
    metric = 'nn';
    % metric = 'cosine_dist';
    n_repeats = 1000;
    n_comp = 25;
    dec_acc = zeros(size(pred_feat, 1), 1);
    % n_distr = 1;
    chance = ones(1, n_comp)./[2:n_comp+1];

    for n_distr = 1:n_comp % takes ~35 min to run when n_comp = 50
        
        dec_acc_im = zeros(size(pred_feat, 1), n_repeats);
        
        for im = 1:size(pred_feat, 1) % for each reconstructed image'
            dist = zeros(1, n_distr+1);
            for rep = 1:n_repeats
                
                % sample images
                % ensure all images including the target are in sample_ids
                sample_set = setdiff(1:size(params, 1), im);
                sample_ids = [randsample(sample_set, n_distr, false) im];
                
                
                sample_ims = params(sample_ids, :);
                
                % test image
                target_im = pred_feat(im, :);
                
                % compute distances of 50 images to target
                if strcmp(metric, 'nn')
                    dist = [];
                    dist = vecnorm((sample_ims - target_im)'); % vecnorm finds norm of each column so do it on transpose
                    
                    % find minimum
                    [comp_dist, comp_idx] = min(dist);
                    
                    if sample_ids(comp_idx) == im
                        dec_acc_im(im, rep) = 1;
                    end
                    
                elseif strcmp(metric, 'cosine_dist')
                    dist = [];
                    dist = sample_ims*target_im'; % project target onto samples
                    
                    % find maximum
                    [comp_dist, comp_idx] = max(dist); % larger proj value = smaller angle between vectors
                    
                    if sample_ids(comp_idx) == im
                        dec_acc_im(im, rep) = 1;
                    end
                end
                
            end
            dec_acc(im, n_distr) = sum(dec_acc_im(im, :))/n_repeats;
            
        end
%         disp(['Finished decoding run with ' num2str(n_distr) ' distractors']);
        
    end
       
    sess_acc(ss, :) = mean(dec_acc, 1);
%     keyboard
    disp(['Finished decoding run for session ' num2str(ss)]);
    toc
end




%% 

patID = {'P76CS', 'P79CS', 'P80CS', 'P84CS', 'P85CS'};
% chance = [1/6 1/6 1/6 1/8 1/8 1/8 1/8 1/8];
p_ID = [1 1 1 2 2 2 3 3 4 4 5]';

f = figure;
hold on
plot(chance, '--k', 'LineWidth', 2)

for i = 1:size(sess_acc, 1)
    
    plot(sess_acc(i, :)', 'LineWidth', 2)
    
end
% legend()
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
title({'Imagined decoding accuracy - axis tuned neurons', [num2str(ndim) ' dims']})
filename = [diskPath filesep 'Recall_Task' filesep ['ImaginedRecon_' num2str(ndim) 'dims']];

print(f, filename, '-dpng', '-r0')

%% 
% 
% testResp = {};
% numcells = 0;
% ctr = 1;
% for ss = 1:length(featureMatrices)
%     
%     fMat = featureMatrices{ss};
%     labels = featureLabels{ss};
%     
%     numcells = numcells + size(fMat, 2);
%     
%     for cellIndex = 1:size(fMat, 2)
%         for ll = 1:length(unique(labels))
%             
%             testResp{ctr, 1}(ll) = mean(fMat(labels == ll, cellIndex))*1e3; 
%             
%         end
%         ctr = ctr+1;
%     end
%     
% end
% 
% 
% 
% %% strctResp Testing. 
% 
% 
% setDiskPaths
% % load([diskPath filesep 'Recall_Task' filesep 'AllITCells_500Stim_Im_SigRamp.mat'])
% % sigRampNames = cell2mat(responses(:, 3));
% 
% sessID = {['Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917'],...
%     ['Recall_Task' filesep 'P76CS' filesep 'Recall_Session_2_20210925'],... % want the rasters so using recall session
%     ['Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_3_20210927'],...
%     ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_1_20220330'],...
%     ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_2_20220403'],...
%     ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_3_20220405'],...
%     ['Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_1_20220728'],...
%     ['Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_2_20220731']};
% 
%     test = load([diskPath filesep sessID{1} filesep 'RecallData_NoFreeRec.mat']);
% 
% 
% sigRampOnly = false;
% strctAllCELLS = {};
% for ss = 4%1:length(sessID)
%     
%     load([diskPath filesep sessID{ss} filesep 'RecallData_NoFreeRec.mat']);
%     strctCELL = struct2cell(strctCells');
%     strctCELL = strctCELL';
%     
%     IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
%     strctCells = strctCells(IT_Cells);
%     
%     RecallData.CRTimeCourse = RecallData.CRTimeCourse(IT_Cells, :);
%     
%     RecallData = rmfield(RecallData, 'CRResponses');
%     
%     strctAllCELLS = [strctAllCELLS; strctCELL(IT_Cells, :)];
%     
%     cellNames = cat(1, strctCells(:).Name);
%     
%     if sigRampOnly
%         if ss ~= 2
%             validNames = ismember(cellNames, sigRampNames);
%             imRasters = imRasters(validNames, :);
%             RecallData.CRTimeCourse = RecallData.CRTimeCourse(validNames, :);
%         end       
%     end
%     
%     
%     CRLength = 5000; %ms
%     for cellIndex = 1:length(strctCells)
%         spikeCount = [];
%         stimRaster = {};
%         for stim = 1:length(RecallData.imageIDs)
%             full_CR_psth = RecallData.CRTimeCourse(cellIndex, 1:3);
%             CR_psth = {full_CR_psth{1, 1}(find(RecallData.CROrder == stim), :), full_CR_psth{1, 2}(find(RecallData.CROrder == stim), :), full_CR_psth{1, 3}};
%             CRStimRaster = CR_psth{1, 1};
%             stimRaster{stim, 1} = CRStimRaster;
%             spikeCount(end+1, :) = mean(mean(CRStimRaster(:, RecallData.offsetTones:RecallData.offsetTones+CRLength)))*1e3;
%         end
%         RecallData.CRResponses{cellIndex, 1} = spikeCount;
%         %     test{cellIndex, 1} = spikeCount;
%     end
%     
%     
%     
% end
% 








