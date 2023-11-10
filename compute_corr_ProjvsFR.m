% Script to compute correlations between projection values and fr
% Currently taking the x coordinates and y coordinates of imagined stimuli from strctCells
% This script is to try corr(value_sta_proj(recalled_stim), fr(recalled_stim)) and see if it is different

%% 

setDiskPaths



load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'params_Alexnet_fc6_500Stimuli.mat']); % will create dataParams = 500x50

% task = 'Object_Screening';
task = 'Recall_Task';

if strcmp(task, 'Recall_Task')
    load([diskPath filesep 'Recall_Task' filesep 'AllITCells_500stim_Im_SigRamp']);
    load([diskPath filesep 'Recall_Task' filesep 'AllITResponses_500stim_Im_SigRamp']);
    lbl = 'Imagination';
    
    if ~exist('ovrlap')
        load([diskPath filesep 'Recall_Task' filesep 'SigRampCellsthatReactivate_alpha0.05.mat'])
    end
    
elseif strcmp(task, 'Object_Screening')
    load([diskPath filesep 'Object_Screening' filesep 'MergedITCells_500stim_Scrn_SigRamp']); % loads cells with sig ramps
    lbl = 'Viewing';
end

ovrlap = ~ovrlap;
% reactivated cells only
if exist('ovrlap', 'var')
    strctCells = strctCells(ovrlap);
    responses = responses(ovrlap, :);
    psths = psths(ovrlap, :);
    strctResp = strctResp(ovrlap);
end

options.task = task;
taskPath = [diskPath filesep task];

taskPath = [taskPath filesep 'forPaper'];
if ~exist(taskPath)
    mkdir([taskPath]);
end

separateSides = false;

if separateSides
    left = 1;
    right = 0;
    strctCELL = struct2cell(strctCells');
    strctCELL = strctCELL';
    
    if left
        ITSide = cellfun(@(x) strcmp(x, 'LFFA'), strctCELL(:, 4));        
    elseif right
        ITSide = cellfun(@(x) strcmp(x, 'RFFA'), strctCELL(:, 4));
    end
    
    strctCells = strctCells(ITSide);
    psths = psths(ITSide, :);
    responses = responses(ITSide, :);
    if exist('strctResp', 'var')
        strctResp = strctResp(ITSide);
    end

end


imageIDs = [1:500];
options.screenType = 'Object';

cc = [];
for cellIndex = 1:length(strctCells)
    
    options.ind_train = imageIDs;
    if strcmp(options.task, 'Recall_Task')
        options.recalledStim = strctCells(cellIndex).recalledStim;
        
        options.ScrnResp = strctResp(cellIndex).ScrnResp;
        options.CRResp = strctResp(cellIndex).CRResp;
        
        options.cellIndex = cellIndex;
    end
    cc(cellIndex, :) = returnCorrelationValues(responses{cellIndex, 1}, params, options);
    
end

mean(cc)
%% cdf
e_cdf = 0;
f = figure;
hold on

if e_cdf
    % % ecdf with bounds - try to make this work
    e1 = ecdf(cc(:, 2), 'Bounds', 'on'); % ortho
    e2 = ecdf(cc(:, 1), 'Bounds', 'on');% pref
    grid on
%     plot(e1, x1, 'LineWidth', 2);
%     plot(e2, x2, 'LineWidth', 2);
    
else
    cd1 = cdfplot(cc(:, 2));
    cd2 = cdfplot(cc(:, 1));
    cd1.LineWidth = 2;
    cd2.LineWidth = 2;

    [c_p, x_p, ~, ~, ~] = cdfcalc(cc(:, 1));
    [c_o, x_o, ~, ~, ~] = cdfcalc(cc(:, 2));
end

[h, p] = kstest2(cc(:, 1), cc(:, 2));


% ylim([0 15])
y_lim = ylim;

if e_cdf
    lgnd = legend({'Ortho axis','', '', 'Preferred axis','', ''});
    filename = [taskPath filesep 'ECDFProjvsFR_STAvsOrtho'];

else
    lgnd = legend({'Ortho axis','Preferred axis'});
    filename = [taskPath filesep 'CDFProjvsFR_STAvsOrtho_' lbl];

end

if exist('ovrlap', 'var')
    filename = [filename '_reactivatedCells'];
end


xlabel('x = Correlation value');
% ylabel('No of neurons');

if strcmp(task, 'Object_Screening')
    xlim([-0.2 1]);
    x_lim = xlim;

    text(0.15*0.80, y_lim(2)*0.88, ['p = ' num2str(p)], 'FontSize', 14,'FontWeight', 'bold');
lgnd.Position = [0.327678571428572,0.646428571428573,0.224107142857143,0.08452380952381];
    title({'Correlation of projection value vs firing rate', 'Preferred and Orthogonal axes', 'Viewing'})
    
elseif strcmp(task, 'Recall_Task')
    x_lim = xlim;
lgnd.Position = [0.170982145837375,0.550000001490118,0.291071422610964,0.110714282734053];

    title({'Correlation of projection value vs firing rate', 'Preferred and Orthogonal axes', 'Imagination'})
    text(x_lim(1)*0.80, y_lim(2)*0.88, ['p = ' num2str(p)], 'FontSize', 14,'FontWeight', 'bold');

end

if separateSides
    if left && ~right
        filename = [filename '_left'];
    elseif right && ~left
        filename = [filename '_right'];
    else
        keyboard
    end
end

set(gca, 'FontSize', 14, 'FontWeight', 'bold');
% print(f, filename, '-dpng', '-r300')
%% histogram of differences between ecdfs of pref and ortho axes
%% This is basically the whole function. Improve this code
clear
setDiskPaths
taskCorr = {};
imageIDs = [1:500];
options.screenType = 'Object';
load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'params_Alexnet_fc6_500Stimuli.mat']); % will create dataParams = 500x50

for t = 1:2
    if t == 1
        task = 'Object_Screening';
    elseif t ==2
        task = 'Recall_Task';
    end
    options.task = task;
    
    if strcmp(task, 'Recall_Task')
        
        load([diskPath filesep 'Recall_Task' filesep 'AllITCells_500stim_Im_SigRamp']);
        load([diskPath filesep 'Recall_Task' filesep 'AllITResponses_500stim_Im_SigRamp']);
        lbl = 'Imagination';
        
    elseif strcmp(task, 'Object_Screening')
        load([diskPath filesep 'Object_Screening' filesep 'MergedITCells_500stim_Scrn_SigRamp']); % loads cells with sig ramps
        lbl = 'Viewing';
    end
    
    
    
    
    cc = [];
    for cellIndex = 1:length(strctCells)
        
        options.ind_train = imageIDs;
        if strcmp(options.task, 'Recall_Task')
            options.recalledStim = strctCells(cellIndex).recalledStim;
            
            options.ScrnResp = strctResp(cellIndex).ScrnResp;
            options.CRResp = strctResp(cellIndex).CRResp;
            
            options.cellIndex = cellIndex;
        end
        cc(cellIndex, :) = returnCorrelationValues(responses{cellIndex, 1}, params, options);
        
    end
    
    [c_p, x_p, ~, ~, ~] = cdfcalc(cc(:, 1));
    [c_o, x_o, ~, ~, ~] = cdfcalc(cc(:, 2));
    
    taskCorr{t} = x_p - x_o;
end
%
% make figure
f = figure; 
hold on
histogram(taskCorr{1}, 'FaceColor', [0.5 0 0.4]) % screening
histogram(taskCorr{2}, 'FaceColor', [0.8500 0 0]) % imagination
xlim([0 1])
y_lim = ylim;
x_lim = xlim;

ylim([0 y_lim(2)*1.1])
y_lim = ylim;
m1 = median(taskCorr{1});
m2 = median(taskCorr{2});



text(m1*0.75, y_lim(2)*0.90, {'median' [ '= ' num2str(m1)]}, 'FontSize', 14,'FontWeight', 'bold');
text(m2*0.65, y_lim(2)*0.67, {'median' [ '= ' num2str(m2)]}, 'FontSize', 14,'FontWeight', 'bold');


lgnd = legend({'Viewing', 'Imagination'});
xlabel('Difference between pref and ortho cdfs');
ylabel('No of neurons');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
filename = [diskPath filesep 'Object_Screening' filesep 'Hist_DiffofCDFs_ViewingandIm'];
% print(f, filename, '-dpng', '-r300')

%% histogram of different correlations
f = figure; 
hold on; 

if strcmp(options.task, 'Recall_Task')
    histogram(cc(:, 2), 'BinEdges', -1:0.1:1) % ortho
    histogram(cc(:, 1), 'BinEdges', -1:0.1:1) % pref
elseif strcmp(options.task, 'Object_Screening')
    histogram(cc(:, 2), 'BinEdges', -0.5:0.05:1) % ortho
    histogram(cc(:, 1), 'BinEdges', -0.5:0.05:1) % pref
    xlim([-0.5 1])
end
% test if distributions are significantly different 
[h, p] = kstest2(cc(:, 1), cc(:, 2));


x_lim = xlim;
% ylim([0 15])
y_lim = ylim;
text(x_lim(1)*0.80, y_lim(2)*0.88, ['p = ' num2str(p)], 'FontSize', 14,'FontWeight', 'bold');
lgnd = legend({'Ortho axis', 'Preferred axis'});
xlabel('Correlation value');
ylabel('No of neurons');
title({'Correlation of projection value vs firing rate', lbl})
set(gca, 'FontSize', 14, 'FontWeight', 'bold');

filename = [taskPath filesep 'CorrProjValvsFR_STAvsOrtho_' lbl];
if exist('ovrlap', 'var')
    filename = [filename '_reactivatedCells'];

end

% print(f, filename, '-dpng', '-r300')


%% Helpers
function [cor] = returnCorrelationValues(fr_raw, params, options)

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
%     touse_fr = fr./std(fr); % already mean subtracted
    touse_fr = fr; % already mean subtracted
elseif strcmp(options.task, 'Recall_Task')
    touse_fr = options.CRResp;
%     im_fr = options.ScrnResp;    
%     if ~isequal(fr_raw(options.recalledStim), im_fr)
%         disp(options.cellIndex);
%     end
    touse_fr = touse_fr - mean(touse_fr); % don't want to subtract the mean from viewing...different conditions
%     touse_fr = touse_fr./std(touse_fr);
end
% PCA
COEFF = pca(para_sub_sta);

pc1 = para_sub_sta * COEFF(:,1); % projections on to principal orthogonal axis


if strcmp(options.task, 'Object_Screening')
    cor(1) = corr(value_sta_prj', touse_fr);
    cor(2) = corr(pc1, touse_fr);
elseif strcmp(options.task, 'Recall_Task')
    cor(1) = corr(value_sta_prj(options.recalledStim)', touse_fr); % note this is mean subtracted but not standardized
    cor(2) = corr(pc1(options.recalledStim), touse_fr); % note this is mean subtracted but not standardized
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

