function outstruct = buildStructWrtTime(obj, timestr, label2find, tfind, deltaT)
% BUILDSTRUCTWRTTIME(obj, timestr, label) Find static/dynamic embryos 
%   Give the times, folders, and time uncertainties of all
%   dynamic samples or fixed samples matching the supplied channel 'label', 
%   depending on whether they match our desired classification as 'dynamic'
%   or 'fixed'.
%
%   This is a helper function for the lookupMap class.
%
% Parameters
% ----------
% obj : the class instance of lookup
% timestr : 'dynamic' or 'static'
% label2find : string, label name (ex 'Eve' or 'Runt')
% 
% Returns 
% -------
% outstruct : struct with fields
%   folders = cell ;
%   names = cell ;
%   embryoIDs = cell ;
%   times = list ;
%   uncs = list ;
%   nTimePoints = list ;
%   tiffpages = list ;
%
% NPMitchell 2020

if nargin < 3
    error("Must supply label (a channel to search for) for class method")
elseif nargin < 5
    % Unless specified, we accept ALL timestamps
    deltaT = Inf ;
    tfind = 0 ;
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
includeTiffpages = false ;

% Now build lists
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
            else
                error('Search for static or dynamic embryos. Which do you want? Pass as argument to buildStructWrtTime()')
            end
            
            if do_add
                if isfinite(deltaT)
                    ts = timestamps{jj} ;
                    for qq = 1:length(ts)
                        if ts(qq) > tfind - deltaT && ts(qq) < tfind + deltaT
                            disp('Adding restricted timestamps...')
                            % Add only certain timestamps, if they match
                            % the condition: tfind-deltaT < t < tfind+deltaT
                            folders{dmyk} = substruct.folders{jj} ;
                            names{dmyk} = substruct.names{jj} ;
                            embryoIDs{dmyk} = substruct.embryoIDs{jj} ;
                            timepts{dmyk} = ts(qq) ;
                            uncs{dmyk} = substruct.uncs{jj}(qq) ;
                            tiffpages{dmyk} = qq ;
                            nTimePoints(dmyk) = substruct.nTimePoints(jj) ;
                            dmyk = dmyk + 1 ;
                            includeTiffpages = true ;
                        end
                    end
                else
                    % Add all timestamps, no matter what they are
                    folders{dmyk} = substruct.folders{jj} ;
                    names{dmyk} = substruct.names{jj} ;
                    embryoIDs{dmyk} = substruct.embryoIDs{jj} ;
                    timepts{dmyk} = timestamps{jj} ;
                    uncs{dmyk} = substruct.uncs{jj} ;
                    nTimePoints(dmyk) = substruct.nTimePoints(jj) ;
                    % tiffpages(dmyk) = 1 ;
                    dmyk = dmyk + 1 ;
                end
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

if includeTiffpages
    outstruct.tiffpages = tiffpages ;
end