% script to call PlotReactivationVis_AlongAxis from


setDiskPaths

% set session list
task = 'Recall_Task';
[sessID, ~, recStim] = Utilities.sessionListAllTasks(task, false, true);
sessID{2} = ['Recall_Task' filesep 'P76CS' filesep 'Recall_Session_2_20210925']; % want the rasters so using recall session

alpha = 0.05;

layermat = 'fc6';
stimDir = '500Stimuli';
load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_Alexnet_' layermat '_' stimDir '.mat']]);

% reacrtive cells
load([diskPath filesep task filesep ['ReactiveITCells_alpha' num2str(alpha) '_500Stim_Im_SigRamp.mat']]);
RC = strctCells;
clearvars psths responses strctCells

RC_names = cat(1, RC(:).Name);
RC_ChanNum = cat(1, RC(:).ChannelNumber);

ITOnly = true;
imageIDs = [1:500]';

plotOptions = struct;
plotOptions.MarkerSize = 4;
plotOptions.Fontsize = 20;
plotOptions.LineWidth = 3;

options.task = task;
options.screenType = 'Object';

for ss = 1:length(sessID)
    tic
    if ss ~=2
        load([sessID{ss} filesep 'RecallData_NoFreeRec.mat'])
        load([sessID{ss} filesep 'PsthandResponses.mat'])
        load([sessID{ss} filesep 'ITResponses.mat'])
        strctCELL = struct2cell(strctCells')';
        
        if ITOnly
            IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
            full_strctCELL = strctCELL;
            strctCELL = strctCELL(IT_Cells, :);
            strctCells = strctCells(IT_Cells);
            EncodingTimeCourse =  RecallData.EncodingTimeCourse(IT_Cells, 1:3);
            CRTimeCourse = RecallData.CRTimeCourse(IT_Cells, 1:3);
            
            psths(~IT_Cells, :) = [];
            responses(~IT_Cells, :) = [];
            
        end
        %%
        resp_ids = cell2mat(cellfun(@(x) ~isnan(x), responses(:, 2), 'UniformOutput', false));
        psths(~resp_ids, :) = [];
        responses(~resp_ids, :) = [];
        strctCells = strctCells(resp_ids);
        EncodingTimeCourse = EncodingTimeCourse(resp_ids, 1:3);
        CRTimeCourse = CRTimeCourse(resp_ids, 1:3);
        % load in stimuli
        stimuli = Utilities.readInFiles([sessID{ss} filesep 'stimuliUsedRecall']);
        
        
        stim_ims = cell(length(stimuli), 1);
        
        for cellIndex = 1:length(strctCells)
            
            
            enc_psth = EncodingTimeCourse(cellIndex, 1:3);
            CR_psth = CRTimeCourse(cellIndex, 1:3);
            EncodingOrder = RecallData.EncodingOrder;
            CROrder = RecallData.CROrder;
            offsetEnc = RecallData.offsetEnc;
            offsetTones = RecallData.offsetTones;
            
            options.ind_train = imageIDs;
            if strcmp(options.task, 'Recall_Task')
                options.recalledStim = recStim{ss};
                
                options.ScrnResp = strctResp(cellIndex).ScrnResp;
                options.CRResp = strctResp(cellIndex).CRResp;
                
                options.cellName = strctResp(cellIndex).Name;
            end
            
            % important bit - figure out order along axis
            axToUse = 'sta';
            orderAlongAx = returnOrderAlongAx(responses{cellIndex, 1}, params, options, axToUse);
%             stimuli = stimuli(orderAlongAx);
            for i = 1:length(stimuli)
                stim_ims{i} = imread([stimuli(i).folder filesep stimuli(i).name]);
                stim_ims{i} = repmat(stim_ims{i}, [1, 1, 3]); 
            end
            f = Utilities.Plotting.PlotReactivationVis_AlongAxis(enc_psth, CR_psth, EncodingOrder, CROrder, offsetEnc, offsetTones, stim_ims, orderAlongAx, plotOptions);
            
            if strcmp(axToUse, 'sta')
                sgt = sgtitle({[num2str(strctCells(cellIndex).Name) ' ' strctCells(cellIndex).brainArea], ['Preferred axis']}); % backslash allows you to print the underscore
            elseif strcmp(axToUse, 'ortho')
                sgt = sgtitle({[num2str(strctCells(cellIndex).Name) ' ' strctCells(cellIndex).brainArea], ['Orthogonal axis']}); % backslash allows you to print the underscore
            end
            
            sgt.FontSize = 20;
            sgt.FontWeight = 'Bold';
            
            pathOut = [diskPath filesep task filesep 'EncodingandCR_comparison' filesep 'AlongAxes'];
            if ~exist(pathOut, 'dir')
                mkdir(pathOut)
            end
            
            if sum(ismember(RC_names, strctCells(cellIndex).Name)) ~= 0 && isequal(strctCells(cellIndex).ChannelNumber, RC_ChanNum(find(ismember(RC_names, strctCells(cellIndex).Name))))
                pathOut = [pathOut filesep 'ReactivatedCells_SigRamp'];
                if ~exist(pathOut, 'dir')
                    mkdir(pathOut)
                end
            end
            filename = [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name) '_Comparison'];
            
            print(f, filename, '-dpng','-r0');
            close all
        end
    end
end

%% helpers
function stimOrd_ax = returnOrderAlongAx(fr_raw, params, options, axis)
% This is a function to


para = params(options.ind_train,:);
amp_dim = sqrt(sum(para.^2)); % finding the norm of each dimension 1xndim
if strcmp(options.screenType, 'Object')
    para = param_normalize_per_dim(para, amp_dim, length(options.ind_train));
elseif strcmp(options.screenType, 'Face')
    para = param_normalize(para, amp_dim, ndim1);
end

ndim = size(para, 2);

fr = fr_raw - mean(fr_raw);
sta=fr'*para;

value_sta_prj = (sta/norm(sta))*para';


para_sub_sta = zeros(size(para));
for k=1:size(para,1);
    param_sta_prj = sta*(para(k,:)*sta')/(sta*sta'); % vector of params pojected onto STA
    para_sub_sta(k,:) = para(k,:) - param_sta_prj; % subtract STA component from param
end

% Note: standardizing makes no difference here
if strcmp(options.task, 'Object_Screening')
    touse_fr = fr; % already mean subtracted
elseif strcmp(options.task, 'Recall_Task')
    touse_fr = options.CRResp;
    touse_fr = touse_fr - mean(touse_fr); % don't want to subtract the mean from viewing...different conditions
    %     touse_fr = touse_fr./std(touse_fr);
end
% PCA
COEFF = pca(para_sub_sta);

pc1 = para_sub_sta * COEFF(:,1); % projections on to principal orthogonal axis


if strcmp(options.task, 'Object_Screening')
    if strcmp(axis, 'sta')
        [~, stimOrd_ax] = sort(value_sta_prj');%, touse_fr);e
    elseif strcmp(axis, 'ortho')
        [~, stimOrd_ax] = sort(pc1, touse_fr);
    end
elseif strcmp(options.task, 'Recall_Task')
    if strcmp(axis, 'sta')
        [~, stimOrd_ax] = sort(value_sta_prj(options.recalledStim)');
    elseif strcmp(axis, 'ortho')
        [~, stimOrd_ax] = sort(pc1(options.recalledStim)); % note this is mean subtracted but not standardized
    end
end
end

function param = param_normalize(param, amp_dim, ndim1)
%% normalize shape/appearance separately while keeping the relative amplitude within shape or appearance dimensions
%% stevens way - in the cell paper
ndim = size(param, 2);
% para = para./repmat(amp_dim, [NIMAGE 1]);

param(:,1:ndim1)=param(:,1:ndim1) / sqrt(sum(amp_dim(1:ndim1).^2)) / sqrt(2);
param(:, ndim1+1:ndim)=param(:, ndim1+1:ndim) / sqrt(sum(amp_dim(ndim1+1:ndim).^2)) / sqrt(2);
end

function param = param_normalize_per_dim(param, amp_dim, NIMAGE)
%% normalize each dimension separately 
%% Liang does this only - May2021
param = param./repmat(amp_dim, [NIMAGE 1]);
end
