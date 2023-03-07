%% set paths and define parameters
%% unique to hpc

dbstop if error
diskPath = '/central/groups/adolphslab/vwadia';
boxPath = '/central/groups/adolphslab/vwadia';

taskCodePath = [boxPath filesep 'recallTaskVarun']; % more events described here
addpath(taskCodePath); setTTLCodes;
addpath([diskPath filesep 'Code' filesep 'SFCpackage' filesep 'helpers']);
addpath([diskPath filesep 'Code' filesep 'SFCpackage' filesep 'SFCieeg']);
addpath(genpath([diskPath filesep 'Code' filesep 'osortTextUI']));
