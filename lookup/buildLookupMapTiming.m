function [map] = buildLookupMapTiming(mutantDir, prepend, postpend, timematfn, timematfn2, save_map) 
%BUILDLOOKUPMAPTIMING Principal function for lookupTable Class object construction
%
% Parameters
% ----------
% mutantDir : str
%   directory path where the staining/label dirs containing experiments, each with pullbacks, are located
% prepend : str
%   prepended string in filename of pullback
% postpend : str
%   postpended string in filename of pullback
% timematfn : str
%   name of the .mat file containing 'matchtime' and 'matchtime_unc'
% timematfn2 : str
%   backup name of a .mat file containing 'matchtime' and 'matchtime_unc' if timematfn is not found
% save_map : bool
%   whether to save the map object to disk
%
% Returns
% -------
% map : containers.Map object
%   The map from fixed sample datasets to time estimate
%


% If prepend and postpend are not defined, give them a default value
if nargin < 7 
    save_map = false ;
    if nargin < 6
        timematfn2 = 'timematch_eve_tmin27_tmax45.mat' ;
        if nargin < 5
            timematfn = 'timematch_EveRunt_tmin27_tmax45.mat' ;
            if nargin < 4
                postpend = '_rot_scaled_view1.tif' ;
                % Similar for prepend
                if nargin < 3
                    prepend = 'Max_Cyl*_2_000001_c' ;
                end
            end
        end
    end
end


%% Build a lookuptable for finding embryos at a given time of a given datatype
map = containers.Map() ;

%% Automatically detect all the labeltype folders in this mutantDir
labelDirs = dir(mutantDir) ;
genoParentDir = labelDirs.folder ;
labelDirs = {labelDirs([labelDirs.isdir]).name} ;
labelDirs = labelDirs(~ismember(labelDirs, {'.', '..', 'figures'})) ;

%% For each labelDir, split into different channels
% Get the channel names from the name of each labelDir
allChannels = {} ;
labelMap = {} ;
% runtMap = zeros(length(labelDirs), 1) ;  % runtMap{i} gives index of channel for labelDirs{i}
for j = 1:length(labelDirs)
    chSetDir = labelDirs{j} ;
    channels = strsplit(chSetDir, '_') ;
    labelMap{j} = channels ;
    for i = 1:length(channels)
        channel = channels{i} ;
        if ~any(strcmp(allChannels, channel))
            allChannels{length(allChannels) + 1} = channel ;
        end
    end
end
% Remove 'Blank' from allChannels
allChannels = allChannels(~ismember(allChannels, 'Blank')) ;

%% Examine each channel in turn 
for chii = 1:length(allChannels)
    channel = allChannels{chii} ;
    
    % Make an output directory for this channel only in the mutantDir
    outFigDir = fullfile(mutantDir, ['figures' filesep channel filesep]) ;

    if ~exist(outFigDir, 'dir')
        mkdir(outFigDir)
    end

    % Select selected Channels: debug
    if any(contains(allChannels, channel))
        disp(['Examining channel ' num2str(chii) ': ' channel])

        % Compile list of filenames
        fileNames = {} ;
        channelDirs = {} ;
        channelTimes = [] ;
        channelTimesUnc = [] ;
        channelsForExperiment = {} ;
        channelNumForExperiment = [] ;
        channelNumsForExperiment = {} ;
        
        for i = 1:length(labelDirs)
            % Note that the runt channel is runtMap(i)
            disp(['Seeking channel in labelDir: ', labelDirs{i}])

            % Check if this labelDir contains the current channel
            if contains(labelDirs{i}, channel) 
                labelDir = fullfile(genoParentDir, labelDirs{i}) ;
                
                % Get all channels in this labelDir
                chSetDir = labelDirs{i} ;
                channels_this_labelDir = strsplit(chSetDir, '_') ;
                

                % This labelDir contains the current channel: Which index?
                chstr = "0" ;
                match = 1 ;
                mapii = labelMap{i} ;
                while strcmp(chstr, "0")
                    if strcmp(mapii{match}, channel)
                        chstr = num2str(match) ;
                    else
                        if match > length(mapii)
                            msg = ['No matching channel despite' ...
                                   ' labelDirNames containing string'] ;
                            error(msg) 
                        else
                            match = match + 1;
                        end
                    end
                end
                disp(['Found matching channel in ' labelDirs{i}])
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Get all possible channel numbers
                if any(contains(channels_this_labelDir, 'Blank'))
                    channelnums = [] ;
                    for qq = 1:length(channels_this_labelDir)
                        if any(contains(allChannels, channels_this_labelDir{qq}))
                            channels_this_labelDir{qq} ;
                            channelnums = [channelnums qq] ;
                        end
                    end
%                     labelDirs{i}
%                     error('ppop')
%                     if strcmp('Slp_Blank_Runt', labelDirs{i})
%                         channelnums
%                         error('here')
%                     end
                else
                    channelnums = 1:length(channels_this_labelDir) ;
                end
                disp(['Found all possible channel numbers in ' labelDirs{i}])
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Assign the filename to be seeking
                filename_pattern = [prepend chstr postpend ] ;    
                % runtname_pattern = [prepend runt_chstr postpend ] ;
                % Note: runt_pattern must uniquely identify the runt file name

                %% Automatically detect all the exptDirs in this labelDir
                exptDirs = dir(fullfile(labelDir, '*')) ;
                exptDirNames = {exptDirs([exptDirs.isdir]).name} ;
                exptParentDir = exptDirs.folder ;
                exptDirs = exptDirNames(~ismember(exptDirNames, {'.', '..', 'figures'})) ;
                % filter out dirs with sorted images
                exptDirs_new = {} ; newkk = 1;
                for qqq = 1:length(exptDirs)
                    if isempty(strfind(exptDirs{qqq}, 'sorted'))
                        exptDirs_new{newkk} = exptDirs{qqq} ;
                        newkk = newkk + 1 ;
                    end
                end
                exptDirs = exptDirs_new ;
                

                %% Iterate through each experiment directory
                for expii = 1:length(exptDirs) 
                    exptDir = fullfile(exptParentDir, exptDirs{expii}) ;
                    disp(['Examining exptDir: ' exptDir])

                    % Runt determines the time
                    try
                        load(fullfile(exptDir, timematfn), 'matchtime', 'matchtime_unc')
                    catch
                        load(fullfile(exptDir, timematfn2), 'matchtime', 'matchtime_unc')
                    end

                    % Add this file and timepoint to the list for this
                    % channel
                    fileNames{length(fileNames) + 1} = filename_pattern ;
                    channelDirs{length(channelDirs) + 1} = exptDir ;
                    channelTimes(length(channelTimes) + 1) = matchtime ;
                    channelTimesUnc(length(channelTimesUnc) + 1) = matchtime_unc ;  
                    channelsForExperiment{length(channelsForExperiment) + 1} = channels_this_labelDir ;
                    channelNumForExperiment(length(channelNumForExperiment) + 1) = match ;
                    channelNumsForExperiment{length(channelNumsForExperiment) + 1} = channelnums ;
                end
            end
        end 
        disp('done building fileNames for this channel')
    end

    % We are done with collating all times and filenames for this channel
    [~, inds] = sort(channelTimes) ;
    channelDirs = channelDirs(inds) ;
    channelTimes = channelTimes(inds) ;
    channelTimesUnc = channelTimesUnc(inds) ;
    channelsForExperiment = channelsForExperiment(inds) ;
    channelNumForExperiment = channelNumForExperiment(inds) ;
    channelNumsForExperiment = channelNumsForExperiment(inds) ;

    for qq = 1:length(channelTimes) 
        % Convert to minutes and round
        time = channelTimes(qq) / 60 ;
        if qq == 1
            substruct.times = [time] ;
            substruct.folders = {channelDirs{qq}} ;
            substruct.uncs = [channelTimesUnc(qq) / 60] ;
            substruct.names = {fileNames{qq}} ;
            % Construct exptname from channelDirs{qq}
            exptname = strsplit(channelDirs{qq}, '/') ;
            substruct.exptnames = {exptname{end}}  ;
            substruct.channels = {channelsForExperiment{qq}} ;
            substruct.channelnums = {channelNumsForExperiment{qq}} ;
            substruct.channelnum = [channelNumForExperiment(qq)] ;
            substruct
        else
            substruct.times(length(substruct.times) + 1) = channelTimes(qq) / 60 ;
            substruct.folders{length(substruct.folders) + 1} = channelDirs{qq} ;
            substruct.uncs(length(substruct.uncs) + 1) = channelTimesUnc(qq) / 60 ;
            substruct.names{length(substruct.names) + 1} = fileNames{qq} ;
            % Construct exptname from channelDirs{qq}
            exptname = strsplit(channelDirs{qq}, '/') ;
            substruct.exptnames{length(substruct.exptnames) + 1} = exptname{end} ;
            substruct.channels{length(substruct.channels) + 1} = channelsForExperiment{qq} ;
            substruct.channelnum(length(substruct.channelnum) + 1) = channelNumForExperiment(qq) ;
            substruct.channelnums{length(substruct.channelnums) + 1} = channelNumsForExperiment{qq} ;
            substruct
        end
    end
    % Add this channel info to the map
    map(channel) = substruct ;
end

if save_map
    % Save the lookup Map
    save_map(fullfile(mutantDir, 'lookuptable_containersMap.mat'), 'map')
end
disp('done building map')


