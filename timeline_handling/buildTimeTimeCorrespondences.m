function [ttc, expts, exptIDs, rebuild_ttc] = ...
    buildTimeTimeCorrespondences(expts, exptIDs, hard, ...
    timelineDir, cpathFnBase, dynamic_embryos, options)
% [ttc, expts, exptIDs] = ...
%     buildTimeTimeCorrespondences(expts, exptIDs, hard, ...
%     timelineDir, extn, cpathFnBase, overwrite, dynamic_embryos)
%
% Parameters
% Returns
% -------
% ttc : cell of cells
%   ttc{ii}{jj} is an Nx2 array where ttc{ii}{jj}(:, 2) is timeline indices
%   for ii and ttc{ii}{jj}(:, 1) is the corresponding placement into
%   timeline jj, as floats.   
% expts : 
% exptIDs :
% options : struct with fields
%   overwrite : bool (default=false)
%       overwrite previous results on disk
%   useOffDiagonalsOnly : bool (default=false)
%       perform all-to-all but only i->j, not j->i, so that we have
%       correspondences
%       1-1 1-2 1-3 .... 1-N
%           2-2 2-3 .... 2-N
%               3-3 .... 3-N
%          ....     ....
%                   .... N-N 
%       If this is false, then we have correspondences
%       1-1 1-2 1-3 .... 1-N
%       2-1 2-2 2-3 .... 2-N
%       3-1 3-2 3-3 .... 3-N
%          ....     ....
%       N-1 N-2 N-3 .... N-N 
%
% NOTE:
% extn specifies the subsampling and method, like:
% extn = [sprintf('_ss%02d', ssfactor) '_' corr_method] ;
% NOTE :
% previously, we defined cpathfn as:
% ijstr = [ '_' eIDi '_' eIDj extn ] ;
% cfn = fullfile(corrOutDir, ['corr' ijstr '.mat']) ;
% disp(['Seeking cfn = ' cfn])
% cpathfn = fullfile(corrPathOutDir, ['cpath' ijstr '.mat']) ;
% Now, cpathfn = sprintf(cpathFnBase, eIDi, eIDj), so
% cpathFnBase = fullfile(corrPathOutDir, ['cpath_%s_%s' extn '.mat'])


% {tj} = {dtj} + frame
% Fit line of best fit to each ridge extraction using time-time-corr cell


% Default options 
overwrite = false ;
useOffDiagonalsOnly = false ;

% Unpack supplied options
if nargin > 6
    if isfield(options, 'overwrite')
        overwrite = options.overwrite ;
    end
    if isfield(options, 'useOffDiagonalsOnly')
        useOffDiagonalsOnly = options.useOffDiagonalsOnly ;
    end
end

% Build ttc: time-time correlation cell
ttcfn = fullfile(timelineDir, 'ttc.mat') ;
if exist(ttcfn, 'file') && ~overwrite
    tmp = load(ttcfn, 'ttc', 'exptIDs', 'hard', 'hard_reference_ID') ;
    ttc = tmp.ttc ;
    
    % Ensure that the exptIDs corresponding to the indices i of i_tau0j are
    % loaded from the mat file itau0jfn
    if ~isfield(tmp, 'exptIDs') || ~isfield(tmp, 'hard')
        % this is bad: the exptIDs were not saved with itau0jfn, so we
        % cannot for certain know that the indices i in i_tau0j correctly
        % index exptIDs. The solution is to erase the file on disk which is
        % missing exptIDs so we can recreate it.
        disp('exptIDs or hard reference index were not saved with itau0jfn!')
        msg = ['We must recreate it!', ...
            'Please erase ', ttcfn, ' from disk now and rerun'];
        error(msg)
    end
    
    % Check that the exptIDs on disk are the same as the ones in RAM
    if ~isequal(exptIDs, tmp.exptIDs)
        % error('handle case here')
        % Rebuild expts in RAM to match the saved one
        for ii = 1:length(exptIDs)
            matchID = find(strcmpi(exptIDs, tmp.exptIDs{ii})) ;
            [filepath,~] = fileparts(dynamic_embryos.folders{matchID}) ;
            expts{ii} = fullfile(filepath, exptIDs{ii}) ;
        end
    end

    exptIDs = tmp.exptIDs ;
    hard = tmp.hard ;
    hard_reference_ID = tmp.hard_reference_ID ;
else
    % Build time-time-correspondence array from scratch
    ttc = cell(length(expts), 1) ;
    for ii = 1:length(expts)
        % create the cell
        ttc{ii} = cell(length(expts), 1) ;
    end
    % Populate the correspondences into cell from mat files
    for ii = 1:length(expts)
        if isnumeric(exptIDs{ii})
            disp(['ref dataset ii = ', num2str(exptIDs{ii})])
        else
            disp(['ref dataset ii = ', exptIDs{ii}])
        end
        for jj = 1:length(expts)
            if isnumeric(exptIDs{jj})
                disp(['current dataset jj = ' num2str(exptIDs{jj})])
            else
                disp(['current dataset jj = ' exptIDs{jj}])
            end
            % Define the correlation (not correspondence) matrix filename
            % This should contain either 'time_correspondences' or
            % 'corrPath' as variables in the .mat file.
            cpathfn = sprintf(cpathFnBase, exptIDs{ii}, exptIDs{jj}) ;
            if exist(cpathfn, 'file')
                try
                    tmp = load(cpathfn, 'time_correspondences');
                    time_correspondences = tmp.time_correspondences ;
                catch
                    tmp = load(cpathfn) ;  
                    time_correspondences = tmp.corrPath(:, [2, 1]) ;
                end
                % Note: time_correspondences(:, 1) is timeline indices for ii
                %   and time_correspondences(:, 2) is corresponding placement into jj
                ttc{ii}{jj} = time_correspondences ;
            else
                disp(['No correspondence for ' num2str(ii) ' to ' num2str(jj)])
            end
        end
    end
    
    % Save the time-time corr cell, along with the exptIDs it indexes
    hard_reference_ID = exptIDs{hard} ;
    save(ttcfn, 'ttc', 'exptIDs', 'hard', 'hard_reference_ID')
end

% Plot the dynamics
close all
colorset = define_colors() ;
for use_offset = [true false]
    offset = 0 ;
    for ii = 1:length(ttc)
        subplot(round(length(ttc)*0.5)+1, 2, ii)
        
        if useOffDiagonalsOnly
            % Do earlier dataset correspondences
            for jj = 1:ii-1
                if ~isempty(ttc{jj}{ii})
                    % displace vertically to match
                    if use_offset
                        i0 = ttc{jj}{ii}(1, 2) ;
                        j0 = ttc{jj}{ii}(1, 1) ;
                        offset = i0 - j0 ;
                    end
                    plot(ttc{jj}{ii}(:, 2), ...
                        ttc{jj}{ii}(:, 1) + offset,...
                        '.-', 'color', colorset(jj, :))
                    hold on;
                end
            end

            % Do this dataset to itself
            for jj = ii 
                if ~isempty(ttc{ii}{jj})
                    plot(ttc{ii}{jj}(:, 1), ttc{ii}{jj}(:, 2), '.', ...
                        'color', colorset(jj, :))
                    hold on;
                end
            end

            % Do later dataset correspondences
            for jj = ii+1:length(ttc)
                if ~isempty(ttc{ii}{jj})
                    % displace vertically to match
                    if use_offset
                        i0 = ttc{ii}{jj}(1, 1) ;
                        j0 = ttc{ii}{jj}(1, 2) ;
                        offset = i0 - j0 ;
                    end
                    plot(ttc{ii}{jj}(:, 1), ...
                        ttc{ii}{jj}(:, 2) + offset,...
                        '.-', 'color', colorset(jj, :))
                    hold on;
                end
            end
        else
            
            % Do all correspondences
            for jj = 1:length(ttc)
                if ~isempty(ttc{ii}{jj})
                    % displace vertically to match
                    if use_offset
                        i0 = ttc{ii}{jj}(:, 1) ;
                        j0 = ttc{ii}{jj}(:, 2) ;
                        offset = mean(i0 - j0) ;
                    end
                    plot(ttc{ii}{jj}(:, 1), ...
                        ttc{ii}{jj}(:, 2) + offset,...
                        '.-', 'color', colorset(jj, :))
                    hold on;
                end
            end            
        end

        % labels
        if isnumeric(exptIDs{ii})
            ylabel(['$\tau(t_{' num2str(exptIDs{ii}) '})$' ], ...
                'Interpreter', 'Latex')
        else
            ylabel(['$\tau(t_{' exptIDs{ii} '})$' ], 'Interpreter', 'Latex')
        end
        if ii == 5 || ii == 6
            xlabel('time, $t_i$', 'Interpreter', 'latex')
        end
        xlim([0, Inf])
        ylim([0, Inf]) 
    end

    % Legend
    subplot(round(length(ttc)*0.5)+1, 2, length(ttc) + 1)
    labels = {} ;
    for ii = 1:length(ttc)
        plot([-1],[-1],'.-', 'color', colorset(ii, :))
        labels{ii} = num2str(ii) ;
        hold on
    end
    legend(labels, 'orientation', 'horizontal')
    xlim([0, 1])
    ylim([0, 1])
    axis off
    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [0 0 16 30]) ;
    set(gca,'fontsize', 12);
    
    if use_offset
        figurefn = fullfile(timelineDir, 'time_correspondences_offset.pdf') ;
    else
        figurefn = fullfile(timelineDir, 'time_correspondences.pdf') ;
    end
    
    saveas(gcf, figurefn)
    
    close all
end
