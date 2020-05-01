function [map] = buildLookupMapTiming(genoDir, prepend, exten, timematfn, save_map) 
%BUILDLOOKUPMAPTIMING Principal function for lookupTable Class object construction
%
% Parameters
% ----------
% genoDir : str
%   directory path where the staining/label dirs containing experiments, each with pullbacks, are located
% prepend : str
%   prepended string in filename of pullback
% exten : str
%   extened string in filename of pullback
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


% If prepend and exten are not defined, give them a default value
if nargin < 6
    save_map = false ;
    if nargin < 5
        timematfn = 'timematch_RuntNanobody_stripe7.mat' ;
        if nargin < 4
            exten = '_rot_scaled_view1.tif' ;
            % Similar for prepend
            if nargin < 3
                prepend = 'Max_Cyl*_2_000001_c' ;
            end
        end
    end
end


%% Build a lookuptable for finding embryos at a given time of a given datatype
map = containers.Map() ;

%% Automatically detect all the labeltype folders in this genoDir
labelDirs = dir(genoDir) ;
genoParentDir = labelDirs.folder ;
labelDirs = {labelDirs([labelDirs.isdir]).name} ;
labelDirs = labelDirs(~ismember(labelDirs, {'.', '..', 'figures'})) ;

%% For each labelDir
% Get the channel names from the name of each labelDir
% allChannels = {} ;
% labelMap = {} ;
% % runtMap = zeros(length(labelDirs), 1) ;  % runtMap{i} gives index of channel for labelDirs{i}
% for j = 1:length(labelDirs)
%     chSetDir = labelDirs{j} ;
%     channels = strsplit(chSetDir, '_') ;
%     labelMap{j} = channels ;
%     for i = 1:length(channels)
%         channel = channels{i} ;
%         if ~any(strcmp(allChannels, channel))
%             allChannels{length(allChannels) + 1} = channel ;
%         end
%     end
% end
% % Remove 'Blank' from allChannels
% allChannels = allChannels(~ismember(allChannels, 'Blank')) ;

%% Examine each fluorescent label in turn 
for ii = 1:length(labelDirs)
    label = labelDirs(ii).name ;
    labelDir = fullfile(genoParentDir, labelDirs{ii}) ;
   
    % Quick filter to ensure that this is a real fluorlabel name
    is_real_label = ~contains(label, 'sorted') ;
    is_real_label = is_real_label && isempty(strfind(exptDirs{qqq}, 'alignment')) ;
    is_real_label = is_real_label && length(label) > 2 ;
    if is_real_label
        disp(['Examining label ' num2str(ii) ': ' label])

        % Compile list of filenames
        fileNames = {} ;
        embryoTimes = [] ;
        embryoTimesUnc = [] ;
        
        % Cycle through all embryos in this labelDir
        embryos = dir(labelDir) ;
        embryoDirs = {labelDirs([labelDirs.isdir]).name} ;
        embryoDirs = labelDirs(~ismember(labelDirs, {'.', '..', 'figures'})) ;

        for ee = 1:length(embryoDirs)
            % Obtain embryo's datestamp
            embryo = embryoDirs{ee} ;
            disp(['Examining embryo ' embryo ' in labelDir: ', labelDirs{ee}])

            % Assign the filename to be seeking
            filename_pattern = [prepend exten ] ;  
            fn = fullfile(labelDir, embryo, filename_pattern) ;
            fnmatch = dir(fn) ;
            if isempty(fnmatch) 
                error(['Could not find data: ' fn])
            else
                filename = fnmatch(1).name;
            end
            
            % Runt determines the time
            try
                load(fullfile(exptDir, timematfn), 'matchtime', 'matchtime_unc')
            catch
                disp(['Could not load timestamp for experiment: ' embryo])
                matchtime = NaN ;
                matchtime_unc = NaN ;
            end

            % Add this file and timepoint to the list for this
            % channel
            fileNames{length(fileNames) + 1} = filename ;
            embryoDirs{length(embryoDirs) + 1} = fullfile(labelDir, embryo) ;
            embryoIDs{length(embryoIDs) + 1} = embryo ;
            embryoTimes(length(embryoTimes) + 1) = matchtime ;
            embryoTimesUnc(length(embryoTimesUnc) + 1) = matchtime_unc ;  
        end 
        disp('done building fileNames for this channel')
    end

    % We are done with collating all times and filenames for this channel
    [~, inds] = sort(embryoTimes) ;
    embryoDirs = embryoDirs(inds) ;
    embryoIDs = embryoIDs(inds) ;
    embryoTimes = embryoTimes(inds) ;
    embryoTimesUnc = embryoTimesUnc(inds) ;

    for qq = 1:length(embryoTimes) 
        % Convert to minutes 
        time = embryoTimes(qq) / 60 ;
        time_unc = embryoTimesUnc(qq) / 60 ;
        if qq == 1
            substruct.times = [time] ;
            substruct.folders = {embryoDirs{qq}} ;
            substruct.uncs = [time_unc] ;
            substruct.names = {fileNames{qq}} ;
            substruct.embryoIDs = {embryoIDs{qq}}  ;
            substruct
        else
            substruct.times(length(substruct.times) + 1) = time ;
            substruct.folders{length(substruct.folders) + 1} = embryoDirs{qq} ;
            substruct.uncs(length(substruct.uncs) + 1) = time_unc ;
            substruct.names{length(substruct.names) + 1} = fileNames{qq} ;
            substruct.embryoIDs{length(substruct.embryoIDs) + 1} = embryoIDs{qq} ;
            substruct
        end
    end
    % Add this channel info to the map
    map(label) = substruct ;
end

if save_map
    % Save the lookup Map
    save_map(fullfile(genoDir, 'lookuptable_containersMap.mat'), 'map')
end
disp('done building map')


