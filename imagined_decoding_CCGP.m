% script to do ccgp for Imagination work
% quick and dirty will have to re-write later.


setDiskPaths

cellType = 5;
% 0 -- All Cells
% 1 -- Resp Cells
% 2 -- SigRampCells
% 3 -- Non-ramp neurons
% 4 -- Non-resp neurons
% 5 -- Reactivated neurons


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
    case 5
%         load([diskPath filesep 'Recall_Task' filesep 'ReactiveITCells_alpha0.05_500Stim_Im_SigRamp_Combo.mat']);
        load([diskPath filesep 'Recall_Task' filesep 'ReactiveITCells_alpha0.05_500Stim_Im_SigRamp.mat']);
        chosenNames = cell2mat(responses(:, 3));
end



sessID = Utilities.sessionListAllTasks('Recall_Task', false);
% sessID{2} = ['Recall_Task' filesep 'P76CS' filesep 'Recall_Session_2_20210925']; % want the rasters so using recall session
%%
tic
strctAllCELLS = {};
for ss = 1:length(sessID)
    if ss ~= 2
        
        % feature matrices are for viewing
        load([diskPath filesep sessID{ss} filesep 'RecallData_NoFreeRec.mat'])
        load([diskPath filesep sessID{ss} filesep 'PsthandResponses.mat'])

        strctCELL = struct2cell(strctCells')';
        
        % carve out IT neurons
        IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
        strctCells = strctCells(IT_Cells);
        responses(~IT_Cells, :) = [];
        psths(~IT_Cells, :) = [];
        imRasters = RecallData.CRTimeCourse(IT_Cells, :);
        
        cellNames = cat(1, strctCells(:).Name);       
        % sanity check
        assert(isequal(cell2mat(responses(:, 3)), cellNames));
        
        
        
        % carving out cells
        if cellType ~= 0            
            validNames = ismember(cellNames, chosenNames);
            
%             disp(sum(validNames));
            
            imRasters = imRasters(validNames, :);
            psths = psths(validNames, :);
            responses = responses(validNames, :); 
            strctCells = strctCells(validNames);
        end
        labels_train = order;
        
        % covert train and test labels to match
        train_files = Utilities.readInFiles([diskPath filesep sessID{ss} filesep 'stimuliUsed'], 'tif');
        train_files = struct2cell(train_files)';
        test_files = Utilities.readInFiles([diskPath filesep sessID{ss} filesep 'stimuliUsedRecall'], 'tif');
        test_files = struct2cell(test_files)';
        train_files = cellfun(@(x) x(1:end-4), train_files(:, 1), 'UniformOutput', false); 
        test_files = cellfun(@(x) x(1:end-4), test_files(:, 1), 'UniformOutput', false); 
        
        matched = find(ismember(train_files, test_files) == 1);
     
        labels_test = RecallData.CROrder;
        
        for id = unique(labels_test)'
            labels_test(labels_test == id) = matched(id);
        end
        
        if cellType ~= 0 
            viewSpikes = cellfun(@(x, y) mean(x(:, floor(y):floor(y)+267), 2), psths(:, 1), responses(:, 2), 'UniformOutput', false);
            viewSpikes = cellfun(@(x) x(ismember(labels_train, matched), :), viewSpikes, 'UniformOutput', false);
        else
            for cellIndex = 1:length(responses)
                if ~isnan(responses{cellIndex, 2})
                    respLat = floor(responses{cellIndex, 2});
                else
                    respLat = 170;
                end
                viewSpikes{cellIndex, 1} = mean(psths{cellIndex, 1}(:, respLat:respLat+267), 2);
                viewSpikes{cellIndex, 1} = viewSpikes{cellIndex, 1}(ismember(labels_train, matched) , :);
            end
            
        end
        labels_train = labels_train(ismember(labels_train, matched));
        viewSpik = [];
        % collect spikes from the entire Im period raster and sum
        % note that when I make strctResp I take mean not sum - April2023
%         imSpikes = cellfun(@(x) sum(x(:, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000), 2), imRasters(:, 1), 'UniformOutput', false);
        imSpikes = cellfun(@(x) mean(x(:, RecallData.offsetTones(1):RecallData.offsetTones(1)+5000), 2), imRasters(:, 1), 'UniformOutput', false);
        imSpik = [];
        
        
        for cellIndex = 1:length(imSpikes)
            imSpik = cat(2, imSpik, imSpikes{cellIndex});
            viewSpik = cat(2, viewSpik, viewSpikes{cellIndex});
        end
        
        
        featMat = struct;
        featMat.Train.Data = viewSpik;
        featMat.Train.Labels = labels_train;
        featMat.Test.Data = imSpik;
        featMat.Test.Labels = labels_test;
        save([diskPath filesep sessID{ss} filesep 'CCGP_featmats_CellType' num2str(cellType)], 'featMat')
        fprintf('Finished saving for session %03d \n', ss)
    end
end
toc


%% load in featmats and decode
setDiskPaths

sessID = Utilities.sessionListAllTasks('Recall_Task', false);
add_PCA = false;
PCA_var_p = 95;
useLDA = false; 
cT_ctr = 5;
for cT = 5%[0 1 2 5]

n_iter = 100;
ctr = 1;
for ss = 6:length(sessID)
    tic
    if ss ~= 2
        load([diskPath filesep sessID{ss} filesep 'CCGP_featmats_CellType' num2str(cT)])
        DataTest = zscore(featMat.Test.Data);
        DataTrain = zscore(featMat.Train.Data);
        
        if ~isempty(DataTrain) && ~isempty(DataTest)
            for i = 1:n_iter
                
                if ~useLDA
                    X = DataTrain;
                    Y = featMat.Train.Labels;
                    
                    t    = templateSVM('Standardize',false,'KernelFunction','linear');
                    Mdl   = fitcecoc(X,Y,'Prior','uniform','Learners',t);
                else
                    
                    if add_PCA
                        
                        [coeff, ~, ~,~, explained] = pca(DataTrain);
                        
                        variance = cumsum(explained);
                        
                        idx_var = find(variance > PCA_var_p,1);
                        %idx_90
                        PCA_DataTrain = DataTrain*coeff;
                        PCA_DataTest = DataTest*coeff;
                        DataTrain = PCA_DataTrain(:,1:idx_var);
                        DataTest = PCA_DataTest(:,1:idx_var);
                    end
                    Mdl = fitcdiscr(DataTrain, featMat.Train.Labels, 'DiscrimType', 'diaglinear');
                end
                predictions = predict(Mdl, DataTest);
                groundTruth = featMat.Test.Labels;
                
                accuracy(i) = sum(predictions==groundTruth)/length(featMat.Test.Labels);
                
            end
            sess_acc(ctr, :) = accuracy;
            ctr = ctr+1;
        end
    end
    fprintf('Finished decoding run for session %03d \n', ss)
    toc
end

sess_acc_byType(:, cT_ctr) = mean(sess_acc, 2);
cT_ctr = cT_ctr + 1;
end

%% plotting
setDiskPaths
load('G:\SUAnalysis\Recall_Task\CCGP_ECOC_Celltypes0125Comboand5.mat')
acc_plot = sess_acc_byType;
acc_plot(:, 4) = []; 
p_chance = [1/6 1/6 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8];
dot_size = 40;

f = figure;
hold on
for i = 1:size(acc_plot, 2)
scatter(1:size(acc_plot, 1), acc_plot(:, i), 'filled')

end
scatter([1:10], p_chance, 2*dot_size, [0 0 0], '_', 'LineWidth', 2)

labels = {'All Cells', 'R Cells', 'Ax Cells', 'Rct Cells'};
lgnd = legend(labels)
lgnd.Position = [0.154166670498394,0.701984119036841,0.21964285331113,0.213095232134774];
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
filename = [diskPath filesep 'Recall_Task' filesep 'CCGP_CellTypes012and5'];
print(f, filename, '-dpng', '-r0')
