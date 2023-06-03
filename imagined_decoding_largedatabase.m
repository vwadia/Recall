% Script to do imagined decoding via linear model
% imaged decoding via linear model


%% Linear model

% compile response matrix for each session
% trials x cells, add a columns of ones
% regress(params, response matrix of n-1 ims)

setDiskPaths
taskPath = 'Recall_Task';

% creates feature matrices and labels
load([diskPath filesep taskPath filesep 'FeatureMatrices_SigRampCells_avgWholePeriod_Im.mat'])
load([diskPath filesep taskPath filesep 'AllITResponses_500Stim_Im_SigRamp.mat'])
load([diskPath filesep taskPath filesep 'AllITCells_500Stim_Im_SigRamp.mat'])

fullnames = cat(1, strctResp(:).Name);

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

% image features
stimDir = '500Stimuli';
load([diskPath filesep 'ObjectSpace' filesep stimDir filesep 'scores_alexnet_comparison.mat']); % will make cell array with fc6 fc6 after relu fc7 fc7 after relu
params_full = scores{1};
[coeff, score, ~, ~, explained, mu] = pca(params_full);

load([diskPath filesep 'Object_Screening' filesep 'resp_AlexNet_fc6BeforeReLu_17856ims']) 
proj_into_500 = (resp_wholeimage - repmat(mu, [size(resp_wholeimage, 1) 1]))*coeff; % coeff = PCs of 500 object space




%%
for ss = 1:length(featureMatrices)
    
    % create resp_mat
    load([diskPath filesep sessID{ss} filesep 'RecallData_NoFreeRec.mat'])
    strctCELL = struct2cell(strctCells');
    strctCELL = strctCELL';
    
    IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
    strctCells = strctCells(IT_Cells);
    names = cat(1, strctCells(:).Name);
    
    num_cells = sum(ismember(fullnames, names));
       
    resp_mat = [];
    resp_mat = cat(2, strctResp(ismember(fullnames, names)).CRResp);
    
    % adding column of ones
    resp_mat = [resp_mat ones(size(resp_mat, 1), 1)];

    all_inds = [1:size(resp_mat, 1)]';
    n_stim = length(all_inds);
    
    ndim = ceil(num_cells*0.75); % need this to be less than the number of images
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
    
    % how to project images???
    proj = proj_into_500(:, 1:ndim); % now the images are projected into the space built by my 500 ims
%     proj = score(setdiff(1:500, recalledStim{ss}), 1:ndim); % now the images are projected into the space built by my 500 ims
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
    
    norm_dist = []; %norm(v_recon - v_original)/norm(v_bestpossrecon - v_original)

    % compute normalized distance
    for im = 1:n_stim
        
        %     if ismember(im, confusing_inds)
        %         keyboard
        %     end
        
        v_orig = params(im, :); % parameters of original image
        
        v_recon = para_recon(im, :);
        
        v_bestPossibleRecon = para_bpr(im, :); % is this right?
        
        
        norm_dist(im) = norm(v_recon - v_orig)/norm(v_bestPossibleRecon - v_orig);
        
    end
    
    clc
    sess_dist{ss} = norm_dist;
    
%     clc
end
%% plot
cols = Utilities.distinguishable_colors(length(sess_dist));
figure; 
hold on
title(['Normalized Distance from ' num2str(size(proj, 1)) ' images']);
for ss = 6%1:length(sess_dist)
    
   plot(sess_dist{ss}, '-o', 'Color', cols(ss, :))
%     keyboard
    
end




