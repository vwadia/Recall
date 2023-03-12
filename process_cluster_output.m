% process cluster outout 
% merges all the output files produced by the various workers on the cluster

% note this contains many anal checks and is very inefficient but I am
% taking *ZERO* chances with my calculations being wrong because variable
% mixes

% vwadia March2023
setDiskPaths

% outputDir = [diskPath filesep 'Recall_Task' filesep 'ppc_log' filesep 'Data' filesep 'HPCOutput'];
% outputDir = [diskPath filesep 'Recall_Task' filesep 'ppc_log' filesep 'Data' filesep 'HPCOutput' filesep 'ppc_LFFACell_LHLFP_sigRamp_SI_Cell'];
% outputDir = [diskPath filesep 'Recall_Task' filesep 'ppc_log' filesep 'Data' filesep 'HPCOutput' filesep 'ppc_LHCell_LFFALFP_sigRamp_SI_Cell'];
% outputDir = [diskPath filesep 'Recall_Task' filesep 'ppc_log' filesep 'Data' filesep 'HPCOutput' filesep 'ppc_RFFACell_RHLFP_sigRamp_SI_Cell'];
% outputDir = [diskPath filesep 'Recall_Task' filesep 'ppc_log' filesep 'Data' filesep 'HPCOutput' filesep 'ppc_RHCell_RFFALFP_sigRamp_SI_Cell'];

inputDirs = [diskPath filesep 'Recall_Task' filesep 'ppc_log' filesep 'Data' filesep 'HPCOutput'];

inDirs = Utilities.readInFiles(inputDirs);
% keeping only desired folders
inDirs = inDirs(~ismember({inDirs.name}, {'combined', 'alreadycombined'}));
inDirs = inDirs(cat(1, inDirs(:).isdir));

for iD = 1:length(inDirs)
    outputDir = [inDirs(iD).folder filesep inDirs(iD).name];
    worker_output = Utilities.readInFiles(outputDir, 'mat');
    ppc_full = {};
    ppc_combined = {};
    tic
    for wrk = 1:length(worker_output)
        
        load([worker_output(wrk).folder filesep worker_output(wrk).name]);
        if wrk == 1
            cellArea = params.cellArea;
            lfpArea = params.lfpArea;
            cellType = params.cells;
            chanType = params.chanType;
            ppc_combined = ppc(:);
            %         ppc_combined_boot = ppc_boot(:);
        end
        nsess = size(ppc, 1);
        
        if ~strcmp(cellArea, params.cellArea) || ~strcmp(lfpArea, params.lfpArea)...
                || ~strcmp(cellType, params.cells) || ~strcmp(chanType, params.chanType)
            error("Params of sessions being combined aren't consistent")
        else
            cellArea = params.cellArea;
            lfpArea = params.lfpArea;
            cellType = params.cells;
            chanType = params.chanType;
        end
        
        % read in ppc and if it is the 4D group then
        ppc_full = [ppc_full ppc(:)];
        ppc_reshaped = ppc(:);
        %     ppc_reshaped_boot = ppc_boot(:);
        if size(ppc_full,2) > 1
            
            md = cell2mat(cellfun(@(x) ndims(x), ppc_full, 'UniformOutput', false));
            %         assert(unique(int(corr(md))) == 1, 'Not all columns are the same - ppc computation is wrong');
            dim = max(md(:)); % number of dimensions there are
            
            toCombine = find(ismember(md(:, 1), max(md)) == 1)';
            
            for row = toCombine
                ppc_combined{row} = cat(dim, ppc_combined{row}, ppc_reshaped{row});
                %             ppc_combined_boot{row} = cat(dim+1, ppc_combined_boot{row}, ppc_reshaped_boot{row});
            end
            
            % average here?
        end
    end
    toc
    
    slashPos = strfind(fname, '/');
    fn = fname(slashPos(end)+1:end);
    dashPos = strfind(fn, '_');
    suffix = fn(dashPos(2)+1:end);
    
    % replace with combined versions
    ppc = reshape(ppc_combined, [nsess 2]);
    % ppc_boot = reshape(ppc_combined_boot, [nsess 2]);
    
    outPath = [diskPath filesep 'Recall_Task' filesep 'ppc_log' filesep 'Data' filesep 'HPCOutput' filesep 'combined'];
    if ~exist(outPath, 'dir')
        mkdir(outPath)
    end
    outfilename = [outPath filesep ['ppc_' suffix '_combined']];
    save(outfilename, 'ppc', 'frq', 'params', 'fname', '-v7.3');
    % save(outfilename, 'ppc', 'ppc_boot', 'frq', 'params', 'fname', '-v7.3');
end
%%

% check if shuffles are actually different 
% nsess = size(ppc, 1);
% for nr = 1:nsess
%     if ~isempty(ppc{nr, 1})
%         if ndims(ppc{nr, 1}) > ndims(ppc{nr, 2})
%             dd = 1;
%         else
%             dd = 2;
%         end
%             
%         for it = 1:size(ppc{nr, dd}, 4)
%             mat = ppc{nr, dd}(:, :, :, it);
%             if it > 1
%                 if isequal(mat, ppc{nr, dd}(:, :, :, it-1))
%                     keyboard
%                 end
%                 disp('nope')
%             end
%         end
%     end
% end

