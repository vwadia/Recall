


setDiskPaths
task = 'Recall_Task';

% scale = 'ppc_linear';
scale = 'ppc_log';

HPC = false; balanced = true;

if HPC && balanced
    datPath = [diskPath filesep task filesep scale filesep 'Data' filesep 'HPCOutput' filesep 'combined'];
    %         fn = 'ppc_LFFACell_LHLFP_sigRamp_ScreeningImagination_Cell_combined';
%         fn = 'ppc_LHCell_LFFALFP_sigRamp_ScreeningImagination_Cell_combined';
    %         fn = 'ppc_RFFACell_RHLFP_sigRamp_ScreeningImagination_Cell_combined';
    fn = 'ppc_RHCell_RFFALFP_sigRamp_ScreeningImagination_Cell_combined';
    
    
    load([datPath filesep fn]);
    
%     for nr = 1:size(ppc, 1)
%         if ~isempty(ppc{nr, 1})
%             if ndims(ppc{nr, 1}) > ndims(ppc{nr, 2})
%                 toAverage = 1;
%             else
%                 toAverage = 2;
%             end
%             
%             ppc{nr, toAverage} = mean(ppc{nr, toAverage}, 4);
%         end
%     end
else
    
    if strcmp(scale, 'ppc_linear')
        frq = linspace(2, 100, 98);
        load([diskPath filesep task filesep scale filesep 'ppc_RITCellRHippLFP_allCells.mat']); noSess2 = 1; Side = 'RFFA';
        
    elseif strcmp(scale, 'ppc_log')
        datPath = [diskPath filesep task filesep scale];

        %     load([diskPath filesep task filesep scale filesep 'ppc_RFFACell_RHLFP_sigRamp_EncodingImagination_cellChans']); Side = 'RFFA';
%         load([diskPath filesep task filesep scale filesep 'ppc_RFFACell_RHLFP_sigRamp_ScreeningImagination_Cell']); Side = 'RFFA';
        
%         load([datPath filesep 'ppc_RFFACell_RHLFP_sigRamp_ScreeningImagination_NoNoise']); Side = 'RFFA';
        load([datPath filesep 'ppc_LFFACell_LHLFP_sigRamp_ScreeningImagination_NoNoise']); Side = 'LFFA';
%         load([datPath filesep 'ppc_RFFACell_RALFP_sigRamp_ScreeningImagination_NoNoise']); Side = 'RFFA';
%         load([datPath filesep 'ppc_LFFACell_LALFP_sigRamp_ScreeningImagination_NoNoise']); Side = 'LFFA';
%         load([datPath filesep 'ppc_RFFACell_RACLFP_sigRamp_ScreeningImagination_NoNoise']); Side = 'RFFA';
%         load([datPath filesep 'ppc_LFFACell_LACLFP_sigRamp_ScreeningImagination_NoNoise']); Side = 'LFFA';
%         load([datPath filesep 'ppc_RFFACell_RSMALFP_sigRamp_ScreeningImagination_NoNoise']); Side = 'RFFA';
%         load([datPath filesep 'ppc_LFFACell_LSMALFP_sigRamp_ScreeningImagination_NoNoise']); Side = 'LFFA';
%         load([datPath filesep 'ppc_RFFACell_ROFLFP_sigRamp_ScreeningImagination_NoNoise']); Side = 'RFFA';
%         load([datPath filesep 'ppc_LFFACell_LOFLFP_sigRamp_ScreeningImagination_NoNoise']); Side = 'LFFA';
        
    end
end

%%
% Computing z-scores with bootstrap distributions
thetaRange = frq >= 3 & frq <= 8; % 3 - 8Hz
alphaRange = frq > 8 & frq <= 12; % 8 - 12Hz
betaRange = frq >= 13 & frq < 30; % 13 - 30Hz
lowGammaRange = frq >= 30 & frq <= 60; % 30 - 60Hz
n_freq = length(frq); 
ppc_z = {};

for cond = 1:size(ppc, 2)
    for session = 1:size(ppc, 1)
        if ~isempty(ppc{session, cond})
            ppc_z{session, cond} = nan(size(ppc{session, cond}));
%             ppc_z{session, cond} = ppc_z{session, cond}(:, :, 1:sum(thetaRange));
            for neuron = 1:size(ppc{session, cond}, 1)
                for channel = 1:size(ppc{session, cond}, 2)
                    
                    %                 ppcIndiv = squeeze(ppc{session, cond}(neuron, channel, thetaRange));
%                     ppcIndiv = ppc{session, cond}(neuron, channel, thetaRange); % for all frequencies
%                     ppcDistRange = squeeze(ppc_boot{session, cond}(neuron, channel, :, thetaRange)); 
                    ppcIndiv = ppc{session, cond}(neuron, channel, :); % for all frequencies
                    ppcDistRange = squeeze(ppc_boot{session, cond}(neuron, channel, :, :));
                    
                    for fr = 1:sum(thetaRange)
                        ppcDist = ~isnan(ppcDistRange(fr, :));
                        ppcDist = ppcDistRange(fr, ppcDist);
                        
                        % fit normal distribution
                        [mu, sig] = normfit(ppcDist);
                        
                        % zscore
                        ppcSub = (ppcIndiv(1, 1, fr) - mu)/sig;
                        
                        sigPPC(fr) = ppcSub;
                        
                        
                        %                     ppc_z{session, cond}(neuron, channel, fr) = ppcSub;
                        
                        %                     sigPPC = ppcIndiv(ppcSub > 2);
                        %                     sigF = frq(thetaRange);
                        %                     sigF = sigF(ppcSub > 2);
                    end
                    if  ~isempty(sigPPC(sigPPC > 2))
                        ppc_z{session, cond}(neuron, channel, :) = ppcIndiv;
                        
                    end
                end
            end
        end
    end
end

% because I'm lazy - for plotting
ppc = ppc_z;
% frq = frq(thetaRange);
%% Cluster based parametric statistics

% note that this will fail in sessions with different numbers of neurons
addpath([diskPath filesep 'Code' filesep 'fieldtrip-20200409']) % Jonathan's version

avSess = true;
sessStats = {};
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
        if avSess
            ppc_C1_all = [ppc_C1_all; ppc_C1]; 
            ppc_C2_all = [ppc_C2_all; ppc_C2];            
        else
            sessStats{session} = Utilities.LFP.clusterStats(ppc_C1, ppc_C2, frqs);
        end
    end
end

if avSess
    
    sessStats = Utilities.LFP.clusterStats(ppc_C1_all, ppc_C2_all, frqs);
    
end

%% p-value of paired test over time

for i = 1:n_freq_C1
    
    [~, pvals(i)] = ttest(ppc_C1_all(:, i), ppc_C2_all(:, i));
    
    
    
end

