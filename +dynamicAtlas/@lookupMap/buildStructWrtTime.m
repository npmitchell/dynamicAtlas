function outstruct = buildStructWrtTime(obj, timestr, label2find)
% BUILDSTRUCTWRTTIME(obj, timestr, label) Find static/dynamic embryos 
%   Give the times, folders, and time uncertainties of all
%   live samples matching the supplied channel 'label'
%
% Parameters
% ----------
% obj : the class instance of lookup
% label2find : string, label name (ex 'Eve' or 'Runt')
%
% NPMitchell 2020

if nargin < 3
    error("Must supply label (a channel to search for) for class method")
end

folders = {}; 
names = {}; 
embryoIDs = {} ;
timepts = {}; 
uncs = {};
nTimePoints = [] ;
dmyk = 1 ;
% get all the labels
labels = obj.map.keys ;
% For each label, add struct.label, folders, uncs
outstruct.times = {} ;
outstruct.folders = {} ;
outstruct.uncs = {} ;
outstruct.nTimePoints = [] ;
for ii = 1:length(labels)
    label = labels{ii} ;
    if strcmp(label, label2find)
        substruct = obj.map(label) ;
        timestamps = substruct.times ;
        for jj = 1:length(timestamps)
            % check if this is part of the stain
            ntp = substruct.nTimePoints(jj) ;
            
            do_add = false ;
            % Choose whether to add based on number of timepoints (ntp)
            if strcmp(timestr, 'static')
                if ntp == 1
                    do_add = true ;
                end
            elseif strcmp(timestr, 'dynamic')
                if ntp > 1
                    do_add = true ;
                end
            end
            
            if do_add
                % Add all timestamps, no matter what they are
                folders{dmyk} = substruct.folders{jj} ;
                names{dmyk} = substruct.names{jj} ;
                embryoIDs{dmyk} = substruct.embryoIDs{jj} ;
                timepts{dmyk} = timestamps{jj} ;
                uncs{dmyk} = substruct.uncs{jj} ;
                nTimePoints(dmyk) = substruct.nTimePoints(jj) ;
                dmyk = dmyk + 1 ;
            end
        end
    end
end

% Add to the output struct
outstruct.folders = folders ;
outstruct.names = names ;
outstruct.embryoIDs = embryoIDs ;
outstruct.times = timepts ;
outstruct.uncs = uncs ;
outstruct.nTimePoints = nTimePoints ;