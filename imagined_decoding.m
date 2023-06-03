% Imagined decoding
% finally getting down to this
% leggo.

% Things I will need:
% feature matrices: cells x trials per session
% labels: the ground truth image identities

% Steps I will need:
% destroy correlations between cells?
% Then decode w/ SVMDecoder (cross validated)
% or with fitcdisc (disciminant analysis)?

% Plot scatter with decoding accuracy per session (maybe across time?)
% First just across all 5s period

% Idea for testing reactivation of cells
% if ttest of baseline vs Im period is significant
    
%% Prep data from each session

% Could compute with all cells, all IT cells, all axis tuned cells
tic
setDiskPaths

cellType = 2;
% 0 -- All Cells
% 1 -- Resp Cells
% 2 -- SigRampCells
% 3 -- Non-ramp neurons
% 4 -- Non-resp neurons


switch cellType
    case 1 % carve out responsive neurons
        load([diskPath filesep 'Recall_Task' filesep 'AllRespITCells_500stim_Im.mat'])
        chosenNames = cell2mat(responses(:, 3));
    case 2 % carve out axis tuned neurons
        load([diskPath filesep 'Recall_Task' filesep 'AllITCells_500Stim_Im_SigRamp.mat'])
        chosenNames = cell2mat(responses(:, 3));
    case 3 % carve out *non* axis tuned neurons (resp and non-resp)
        load([diskPath filesep 'Recall_Task' filesep 'AllITCells_500Stim_Im_SigRamp.mat'])
        chosenNames = cell2mat(responses(:, 3));
    case 4 % carve out *non* responsive neurons
        load([diskPath filesep 'Recall_Task' filesep 'AllRespITCells_500stim_Im.mat'])
        chosenNames = cell2mat(responses(:, 3));
end



sessID = Utilities.sessionListAllTasks('Recall_Task', false);
sessID{2} = ['Recall_Task' filesep 'P76CS' filesep 'Recall_Session_2_20210925']; % want the rasters so using recall session


strctAllCELLS = {};
for ss = 1:length(sessID)
    
    load([diskPath filesep sessID{ss} filesep 'RecallData_NoFreeRec.mat'])
    if ss ~= 2
        load([diskPath filesep sessID{ss} filesep 'strctCells'])
    end
    strctCELL = struct2cell(strctCells');
    strctCELL = strctCELL';
    
    IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
    strctCells = strctCells(IT_Cells); 

    cellNames = cat(1, strctCells(:).Name);
    imRasters = RecallData.CRTimeCourse(IT_Cells, :);

    if cellType ~= 0
        
        if cellType == 1 || cellType == 2
            if ss ~= 2
                validNames = ismember(cellNames, chosenNames);
                imRasters = imRasters(validNames, :);
            end
        elseif cellType == 3 || cellType == 4
            if ss == 2
                validNames == zeros(length(cellNames), 1);
            else
                validNames = ~ismember(cellNames, chosenNames);
                imRasters = imRasters(validNames, :);
            end
        end
        
    end
    labels = RecallData.CROrder;
    
    % collect spikes from the entire Im period raster and sum
    % note that when I make strctResp I take mean not sum - April2023
    imSpikes = cellfun(@(x) sum(x(:, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000), 2), imRasters(:, 1), 'UniformOutput', false);
%     imSpikes = cellfun(@(x) mean(x(:, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000), 2), imRasters(:, 1), 'UniformOutput', false);
    imSpik = [];
    
    
    for cellIndex = 1:length(imSpikes)
        imSpik = cat(2, imSpik, imSpikes{cellIndex});
    end
    
    featureMatrices{ss} = imSpik;
    featureLabels{ss} = labels;
    
end
toc
%% decode via SVM

for sess = 1:length(sessID)
    tic
    if ~isempty(featureMatrices{sess})
        % choose a session
        X = featureMatrices{sess};
        t = featureLabels{sess};
        
        X = zscore(X);
        destroyCorr = true;
        
        % destroy correlations across cells
        if destroyCorr
            unique_labels = unique(t);
            X_new         = [];
            for i=1:size(X,2) % for all cells
                temp = [];
                for j=1:size(unique_labels)
                    idx = Utilities.Shuffle(find(ismember(t,unique_labels(j)))); % finds all the scene n labels and shuffles them
                    temp = cat(1,temp,X(idx,i));  % cat(dim, a, b) so dim = 1 = rows. cat(1, ) is the same as vertcat
                end
                X_new(:,i) = temp;
            end
            t = sort(t);
        else
            X_new = X;
        end
        
        % run Decoder
        % X_new = X;
        nr_iterations = 100;
        
        accuracy = [];
        for i=1:nr_iterations % need to run a few times to get stable estimates
            %         idx_train     = Utilities.Decode.balance_groups(t);
            %         X_train       = X_new(cell2mat(idx_train),:);
            %         t_train      = t(cell2mat(idx_train));
            
            % my data is trial balanced so no need to balance groups first?
            X_train = X;
            t_train = t;
            
            [predictions,groundtruth, Mdl,~,~] = Utilities.Decode.SVMdecoder(X_train,t_train); % why does Mdl.Trained have 10 values? = Per k-fold!
            accuracy(i) = sum(predictions==groundtruth)/length(t_train);
            %         fprintf('Finished decoding run %03d \n', i) % looks nicest
        end
        
        fprintf('Finished decoding runs for session %03d \n', sess)
        
        acc(sess, :) = accuracy;
    end
    toc
end


switch cellType
    case 0 % all cells
        save([diskPath filesep 'Recall_Task' filesep 'FeatureMatrices_AllITCells_avgWholePeriod_Im'], 'featureMatrices', 'featureLabels')
        save([diskPath filesep 'Recall_Task' filesep 'ImDecAcc_AllITCells_' num2str(nr_iterations) '_SVM_6or8foldCrossval_PerSession'], 'acc')
    case 1 % resp cells
        save([diskPath filesep 'Recall_Task' filesep 'FeatureMatrices_RespCells_avgWholePeriod_Im'], 'featureMatrices', 'featureLabels')
        save([diskPath filesep 'Recall_Task' filesep 'ImDecAcc_RespCells_' num2str(nr_iterations) '_SVM_6or8foldCrossval_PerSession'], 'acc')
    case 2 % sigramp
        save([diskPath filesep 'Recall_Task' filesep 'FeatureMatrices_SigRampCells_avgWholePeriod_Im'], 'featureMatrices', 'featureLabels')
        save([diskPath filesep 'Recall_Task' filesep 'ImDecAcc_SigRampCells_' num2str(nr_iterations) '_SVM_6or8foldCrossval_PerSession'], 'acc')
    case 3 % non-sigramp
        save([diskPath filesep 'Recall_Task' filesep 'FeatureMatrices_NONSigRampCells_avgWholePeriod_Im'], 'featureMatrices', 'featureLabels')
        save([diskPath filesep 'Recall_Task' filesep 'ImDecAcc_NONSigRampCells_' num2str(nr_iterations) '_SVM_6or8foldCrossval_PerSession'], 'acc')
    case 4 % non-resp
        save([diskPath filesep 'Recall_Task' filesep 'FeatureMatrices_NONRespCells_avgWholePeriod_Im'], 'featureMatrices', 'featureLabels')
        save([diskPath filesep 'Recall_Task' filesep 'ImDecAcc_NONRespCells_' num2str(nr_iterations) '_SVM_6or8foldCrossval_PerSession'], 'acc')
end



%% decode via LDA
setDiskPaths
load([diskPath filesep 'Recall_Task' filesep 'FeatureMatrices_SigRampCells_avgWholePeriod_Im'])

for sess = 1:length(sessID)
    
    % choose a session
    X = featureMatrices{sess};
    t = featureLabels{sess};
    
    parameters = struct;
    parameters.NbrRep = 10;
    parameters.PCA_var_p = 95;
    [acc_train{sess, 1}, acc_test{sess, 1}] = Utilities.Decode.classification_PCA_LDA(X, t, parameters);
    
    
end

% 
% for sess = 1:length(sessID)
%     acc_LDA_PCA(sess, 1) = mean(acc_train{sess, 1});
%     acc_LDA_PCA(sess, 2) = mean(acc_test{sess, 1});
% end
% save([diskPath filesep 'Recall_Task' filesep 'ImDecAcc_' num2str(parameters.PCA_var_p) 'Var_PCA_LDA_6or8foldCrossval_PerSession'], 'acc_LDA_PCA')
% 


%% Per patient plot

sess_acc = mean(acc, 2);
patID = {'P76CS', 'P79CS', 'P80CS', 'P84CS', 'P85CS'};
chance = [1/6 1/6 1/6 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8];
p_ID = [1 1 1 2 2 2 3 3 4 4 5]';
for i = unique(p_ID)'
    p_chance(i) = mean(chance(find(p_ID == i)));
end
% cols = ([0 0.4470 0.7410; 0.4660 0.6740 0.1880; 0.8500 0.3250 0.0980; ]);
cols = Utilities.distinguishable_colors(length(patID));

pat_acc = [mean(sess_acc(1:3)) mean(sess_acc(4:6)) mean(sess_acc(7:8)) mean(sess_acc(9:10)) mean(sess_acc(11))];

dot_size = 40;
f = figure; 
hold on
gscatter(p_ID,  sess_acc, p_ID, cols)
% scatter([1 2 3], pat_acc, dot_size, [0 0 0], '_', 'LineWidth', 2)
scatter([1 2 3 4 5], p_chance, 2*dot_size, [0 0 0], '_', 'LineWidth', 2)
ylim([0 1])
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Decoding accuracy')
xlabel('Patient ID')
lgnd = legend(patID);
lgnd.Position = [0.687440477627374,0.577807542964343,0.187499997126205,0.264285706835134];

xticks([1 2 3 4 5])
% xlim([0 4])

switch cellType
    case 0
        filename = [diskPath filesep 'Recall_Task' filesep 'ImaginedDecoding_AllITCells_avgAcrossWholePeriod'];
        title({'Per session imagined decoding accuracy', 'All IT Neurons'})
    case 1
        filename = [diskPath filesep 'Recall_Task' filesep 'ImaginedDecoding_RespCells_avgAcrossWholePeriod'];
        title({'Per session imagined decoding accuracy', 'Visually Responsive Neurons'})
    case 2
        filename = [diskPath filesep 'Recall_Task' filesep 'ImaginedDecoding_SigRampCells_avgAcrossWholePeriod'];
        title({'Per session imagined decoding accuracy', 'Axis tuned Neurons'})
    case 3
        filename = [diskPath filesep 'Recall_Task' filesep 'ImaginedDecoding_NONSigRampCells_avgAcrossWholePeriod'];
        title({'Per session imagined decoding accuracy', 'Non Axis tuned Neurons'})
    case 4
        filename = [diskPath filesep 'Recall_Task' filesep 'ImaginedDecoding_NONRespCells_avgAcrossWholePeriod'];
        title({'Per session imagined decoding accuracy', 'Non Responsive Neurons'})
end

% print(f, filename, '-dpng', '-r0')


%% Per session plot
% accur = [];
% chance = [];
% for ac = 1:length(acc)
%     
%     accur = cat(1, accur, acc{ac});
%     
%     chance = cat(1, chance, 1/length(unique(featureLabels{ac})));
% 
%     
% end
% dot_size = 40;
% sess_acc = mean(accur, 2);
% figure; 
% hold on
% scatter(1:length(acc), sess_acc, dot_size, [0.6350 0.0780 0.1840], 'filled')
% scatter(1:length(acc), chance, dot_size ,[0 0 0], '_', 'LineWidth', 2)
% ylim([0 1])
% xlim([0 9])
