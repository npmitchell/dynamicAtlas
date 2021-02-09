function [map] = buildLookupMap(lum, options) 
%BUILDLOOKUPMAPTIMING Principal function for lookupTable Class object construction
%
% Parameters
% ----------
% lum : lookupMap class instance with properties
%     genoDir : str
%       directory path where the staining/label dirs containing experiments, each with pullbacks, are located
%     prepend : str
%       prepended string in filename of pullback
%     exten : str
%       extened string in filename of pullback
%     timematfn : str or cell array of strings
%       If string, name of the .mat file or .txt file containing 
%       'matchtime' and 'matchtime_unc'.
%       If cell array, list of names of .mat files or .txt files in rank
%       order to use for time matching
% optionns : struct with fields
%   labels : cell array 
%       the labels to include, if not all labels of a genotype
%   save_map : bool
%       whether to save the map object to disk
%
% Returns
% -------
% map : containers.Map object
%   The map from fixed sample datasets to time estimate
%

% Unpack the supplied lookupMap (lum)
genoDir = lum.genoDir ;
prepend = lum.prepend ;
exten = lum.exten ;
timematfn = lum.timerfn ; % note: could be string for filename match OR 
                          % cell array of strings in rank preference order

% If prepend and exten are not defined, give them a default value
if nargin < 2
    options = struct() ;
end

if isfield(options, 'save_map')
    save_map = options.save_map ;
else
    save_map = false ;
end
if isfield(options, 'labels')
    labels = options.labels ;
else
    labels = {} ;
end

%% Build a lookuptable for finding embryos at a given time of a given datatype
map = containers.Map() ;

%% Automatically detect all the labeltype folders in this genoDir
if isempty(labels)
    disp(['Seeking labels in genotype directory: ' genoDir])
    labelDirs = dir(genoDir) ;
    genoParentDir = labelDirs.folder ;
    labelDirs = {labelDirs([labelDirs.isdir]).name} ;
    labelDirs = labelDirs(~ismember(labelDirs, {'.', '..', 'figures'})) ;
else
    disp(['Seeking only specified labels in genotype directory: ' genoDir])
    labelDirs = dir(genoDir) ;
    genoParentDir = labelDirs.folder ;
    labelDirs = {labelDirs([labelDirs.isdir]).name} ;
    labelDirs = labelDirs(ismember(labelDirs, labels)) ;
end

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
    label = labelDirs{ii} ;
    labelDir = fullfile(genoParentDir, labelDirs{ii}) ;
   
    % Quick filter to ensure that this is a real fluorlabel name
    is_real_label = ~contains(label, 'sorted') ;
    is_real_label = is_real_label && isempty(strfind(labelDirs{ii}, 'alignment')) ;
    is_real_label = is_real_label && length(label) > 2 ;
    if is_real_label
        disp(['Examining label ' num2str(ii) ': ' label])

        % Compile list of filenames
        fileNames = {} ;
        embryoDirs = {} ;
        embryoIDs = {}  ;
        embryoTimes = {} ;
        embryoTimesUnc = {} ;
        nTimePoints = [] ;
        
        % Cycle through all embryos in this labelDir
        embryos = dir(labelDir) ;
        embryos = {embryos([embryos.isdir]).name} ;
        embryos = embryos(~ismember(embryos, {'.', '..', 'figures'})) ;

        for ee = 1:length(embryos)
            % Obtain embryo's datestamp
            embryo = embryos{ee} ;
            disp(['Examining embryo ', embryos{ee}])
            
            % Assign the filename to be seeking
            filename_pattern = [prepend exten ] ;  
            fn = fullfile(labelDir, embryo, filename_pattern) ;
            fnmatch = dir(fn) ;
            if isempty(fnmatch) 
                error(['Could not find data: ' fn])
            else
                filename = fnmatch(1).name;
            end
            
            % Determine number of timepoints  
            ntps = length(imfinfo(fullfile(fnmatch(1).folder, ...
                fnmatch(1).name))) ;
            
            %% Determine the timestamp(s)
            % First convert to cell if not already cell. Iterate over each
            % entry as rank preference for timestamp type to load
            if isa(timematfn, 'char')
                timematfn = {timematfn} ;
            end
            
            if isa(timematfn, 'cell')
                for qq = 1:length(timematfn)
                    tmatMatches = dir(fullfile(fnmatch(1).folder, timematfn{qq})) ;
                    if ~isempty(tmatMatches)
                        matmatch = tmatMatches(1).name ;
                        disp(['loading time from ' matmatch])
                        if strcmpi(matmatch(end-2:end), 'mat')
                            load(fullfile(fnmatch(1).folder, matmatch), 'matchtime_minutes', 'matchtime_unc_minutes')
                        elseif strcmpi(matmatch(end-2:end), 'txt')
                            % Interpret first column as timestamps, second as
                            % time uncertainties
                            matchtime_min_unc = dlmread(fullfile(fnmatch(1).folder, matmatch), ',') ;
                            try
                                assert(all(size(matchtime_min_unc) == [ntps, 2]))
                                matchtime_minutes = matchtime_min_unc(:, 1) ;
                                matchtime_unc_minutes = matchtime_min_unc(:, 2) ;
                            catch
                                error(['Error with file: If timestamps stored in txt file, must be two columns with values, uncertainties: ' matmatch])
                            end
                        else
                            load(fullfile(fnmatch(1).folder, matmatch), 'matchtime_minutes', 'matchtime_unc_minutes')
                        end
                    else
                        disp(['Could not load timestamp for experiment: ' embryo])
                        matchtime_minutes = NaN * ones(1, ntps) ;
                        matchtime_unc_minutes = NaN * ones(1, ntps) ;
                    end
                end                
            else
                error('Could not interpret timematfn=lum.timerfn --> should be string or cell of strings in rank preference order')
            end

            % Add this file and timepoint to the list for this
            % channel
            fileNames{length(fileNames) + 1} = filename ;
            embryoDirs{length(embryoDirs) + 1} = fullfile(labelDir, embryo) ;
            embryoIDs{length(embryoIDs) + 1} = embryo ;
            embryoTimes{length(embryoTimes) + 1} = matchtime_minutes ;
            embryoTimesUnc{length(embryoTimesUnc) + 1} = matchtime_unc_minutes ; 
            nTimePoints(length(nTimePoints) + 1) = ntps ; 
        end 
        disp('done building fileNames for this channel')
    end

    % We are done with collating all times and filenames for this channel
    % Sort them by the minimum time in each Tiff page / stack
    min_etime = zeros(size(embryoTimes)) ;
    for qq = 1:length(min_etime)
        min_etime(qq) = min(embryoTimes{qq}) ;
    end
    
    [~, inds] = sort(min_etime) ;
    embryoDirs = embryoDirs(inds) ;
    embryoIDs = embryoIDs(inds) ;
    embryoTimes = embryoTimes(inds) ;
    embryoTimesUnc = embryoTimesUnc(inds) ;
    nTimePoints = nTimePoints(inds) ;

    for qq = 1:length(embryoTimes) 
        % Convert to minutes 
        time = embryoTimes{qq} ;
        time_unc = embryoTimesUnc{qq} ;
        if qq == 1
            substruct.times = {time} ;
            substruct.folders = {embryoDirs{qq}} ;
            substruct.uncs = {time_unc} ;
            substruct.names = {fileNames{qq}} ;
            substruct.embryoIDs = {embryoIDs{qq}}  ;
            substruct.nTimePoints = [nTimePoints(qq)] ;
        else
            substruct.folders{length(substruct.folders) + 1} = embryoDirs{qq} ;
            substruct.times{length(substruct.times) + 1} = time ;
            substruct.uncs{length(substruct.uncs) + 1} = time_unc ;
            substruct.names{length(substruct.names) + 1} = fileNames{qq} ;
            substruct.embryoIDs{length(substruct.embryoIDs) + 1} = embryoIDs{qq} ;
            substruct.nTimePoints(length(substruct.nTimePoints) + 1) = nTimePoints(qq) ;
        end
    end
    % Add this channel info to the map
    disp(substruct)
    map(label) = substruct ;
end

if save_map
    % Save the lookup Map
    save_map(fullfile(genoDir, 'lookuptable_containersMap.mat'), 'map')
end
disp('done building map')



