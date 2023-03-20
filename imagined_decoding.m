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

setDiskPaths

sessID = {['Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917'],...
    ['Recall_Task' filesep 'P76CS' filesep 'Recall_Session_2_20210925'],... % want the rasters so using recall session
    ['Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_3_20210927'],...
    ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_1_20220330'],...
    ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_2_20220403'],...
    ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_3_20220405'],...
    ['Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_1_20220728'],...
    ['Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_2_20220731']};



for ss = 1:length(sessID)
    
    load([diskPath filesep sessID{ss} filesep 'RecallData_NoFreeRec.mat'])
    strctCELL = struct2cell(strctCells');
    strctCELL = strctCELL';
    
    IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
    strctCells = strctCells(IT_Cells);

    imRasters = RecallData.CRTimeCourse(IT_Cells, :);
    
    labels = RecallData.CROrder;
    
    % collect spikes from the entire Im period raster and sum
    imSpikes = cellfun(@(x) sum(x(:, RecallData.offsetTones(1):RecallData.offsetTones(2)), 2), imRasters(:, 1), 'UniformOutput', false);
    imSpik = [];
    
    for cellIndex = 1:length(strctCells)
        imSpik = cat(2, imSpik, imSpikes{cellIndex});
    end
    
    featureMatrices{ss} = imSpik;
    featureLabels{ss} = labels;
    
    
end

%% decode via SVM

for sess = 1:length(sessID)
    
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

        % my data is trial balanced so no need to balance groups first
        X_train = X;
        t_train = t;
        
        [predictions,groundtruth, Mdl,~,~] = Utilities.Decode.SVMdecoder(X_train,t_train); % why does Mdl.Trained have 10 values? = Per k-fold!
        accuracy(i) = sum(predictions==groundtruth)/length(t_train);
%         fprintf('Finished decoding run %03d \n', i) % looks nicest
    end
    
    fprintf('Finished decoding runs for session %03d \n', sess) 
  
    acc(sess, :) = accuracy;
    
end

save([diskPath filesep 'Recall_Task' filesep 'FeatureMatrices_avgWholePeriod_Im'], 'featureMatrices', 'featureLabels')
save([diskPath filesep 'Recall_Task' filesep 'ImDecAcc_' num2str(nr_iterations) '_SVM_6or8foldCrossval_PerSession'], 'acc')

%% decode via LDA

for sess = 1:length(sessID)
    
    % choose a session
    X = featureMatrices{sess};
    t = featureLabels{sess};
    
    parameters = struct;
    parameters.NbrRep = 10;
    parameters.PCA_var_p = 95;
    [acc_train{sess, 1}, acc_test{sess, 1}] = Utilities.Decode.classification_PCA_LDA(X, t, parameters);
    
    
end


for sess = 1:length(sessID)
    acc_LDA_PCA(sess, 1) = mean(acc_train{sess, 1});
    acc_LDA_PCA(sess, 2) = mean(acc_test{sess, 1});
end
save([diskPath filesep 'Recall_Task' filesep 'ImDecAcc_' num2str(parameters.PCA_var_p) 'Var_PCA_LDA_6or8foldCrossval_PerSession'], 'acc_LDA_PCA')



%% Per patient plot

patID = {'P76CS', 'P79CS', 'P80CS'};

p_ID = [1 1 1 2 2 2 3 3]';
p_chance = [mean(chance(find(p_ID == 1))) mean(chance(find(p_ID == 2))) mean(chance(find(p_ID == 3)))];
cols = ([0 0.4470 0.7410; 0.4660 0.6740 0.1880; 0.8500 0.3250 0.0980]);

pat_acc = [mean(sess_acc(1:3)) mean(sess_acc(4:6)) mean(sess_acc(7:8))];
dot_size = 40;
f = figure; 
hold on
gscatter(p_ID,  sess_acc, p_ID, cols)
% scatter([1 2 3], pat_acc, dot_size, [0 0 0], '_', 'LineWidth', 2)
scatter([1 2 3], p_chance, 2*dot_size, [0 0 0], '_', 'LineWidth', 2)
ylim([0 1])
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Decoding accuracy')
xlabel('Patient ID')
lgnd = legend(patID);
lgnd.Position = [0.682083334770231,0.74566468433137,0.187499997126205,0.161904757434414];
xticks([1 2 3])
title('Per session imagined decoding accuracy')
% xlim([0 4])
filename = [diskPath filesep 'Recall_Task' filesep 'ImaginedDecoding_avgAcrossWholePeriod'];
print(f, filename, '-dpng', '-r0')


%% Per session plot
accur = [];
chance = [];
for ac = 1:length(acc)
    
    accur = cat(1, accur, acc{ac});
    
    chance = cat(1, chance, 1/length(unique(featureLabels{ac})));

    
end
dot_size = 40;
sess_acc = mean(accur, 2);
figure; 
hold on
scatter(1:length(acc), sess_acc, dot_size, [0.6350 0.0780 0.1840], 'filled')
scatter(1:length(acc), chance, dot_size ,[0 0 0], '_', 'LineWidth', 2)
ylim([0 1])
xlim([0 9])
