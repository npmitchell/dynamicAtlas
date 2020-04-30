

clear all
mutantDir = '../Tlrm9/' ;
postpend = '_rot_scaled_view1.tif' ;
prepend = 'Max_Cyl*_2_000001_c' ;
timematfn = 'timematch_EveRunt_tmin27_tmax45.mat' ;
timematfn2 = 'timematch_eve_tmin27_tmax45.mat' ;
save_map = true ;
timematfn = 'timematch_EveRunt_tmin27_tmax45.mat' ;
timematfn2 = 'timematch_eve_tmin27_tmax45.mat' ;

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

% debug:
% allChannels = {'Eve'}

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

        % Compile list of filenames, according to filename_pattern
        fileNames = {} ;
        channelDirs = {} ;
        channelTimes = [] ;
        channelTimesUnc = [] ;
        for i = 1:length(labelDirs)
            % Note that the runt channel is runtMap(i)
            disp(['Seeking channel in labelDir: ', labelDirs{i}])

            % Check if this labelDir contains the current channel
            if contains(labelDirs{i}, channel) 
                labelDir = fullfile(genoParentDir, labelDirs{i}) ;

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
                    channelDirs{length(channelDirs) + 1} = exptDir ;
                    channelTimes(length(channelTimes) + 1) = matchtime ;
                    channelTimesUnc(length(channelTimes) + 1) = matchtime_unc ;  
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

    for qq = 1:length(channelTimes) 
        % Convert to minutes and round
        time = channelTimes(qq) / 60 ;
        if qq == 1
            substruct.times = [time] ;
            substruct.folders = {channelDirs{qq}} ;
            substruct.uncs = [channelTimesUnc(qq) / 60] ;
        else
            substruct.times(length(substruct.times) + 1) = channelTimes(qq) / 60 ;
            substruct.folders{length(substruct.folders) + 1} = channelDirs{qq} ;
            substruct.uncs(length(substruct.uncs) + 1) = channelTimesUnc(qq) / 60 ;
        end
    end
    map(channel) = substruct ;
end

if save_map
    % Save the lookup Map
    save(fullfile(mutantDir, 'lookuptable_containersMap_struct.mat'), 'map')
end
disp('Done with map containing struct.')


