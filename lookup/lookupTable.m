classdef lookupTable
    %LOOKUPTABLE Look up timing or channels from a containers.Map object
    %   Class for a lookuptable for finding times and channels
    %
    % Example
    % -------
    % a = lookupTable
    % a = a.buildLookup('../WT') 
    % a.findTime(40)
    % a.findTime(40, 3)
    % a.findLabel('Eve')
    % a.findLabel('Runt')
    % a.findLabelTime('Eve', 40)
    % a.findLabelTime('Eve', 40, 3)
    %
    % testing
    % -------
    % keys = a.map('Eve').keys
    % for i = 1:length(keys)
    %   q = a.map('Eve') ;
    %   q(keys{i})
    % end
    % 
    % NPMitchell 2019
    
    properties
        
        map = containers.Map
        
    end
    
    
    methods
        function obj = buildLookup(obj, mutantDir, prepend, postpend)
            %BUILDLOOKUP Construct the lookup table
            if nargin < 2
                error('Must supply a directory to use (mutantDir) for class method buildLookup')
            elseif nargin < 3
                prepend = 'Max_Cyl*_2_000001_c' ;       % string before channel index
            end
            if nargin < 4
                postpend = '_rot_scaled_view1.tif' ;    % string after channel index
            end
            timematfn = 'timematch_EveRunt_tmin27_tmax45.mat' ;
            timematfn2 = 'timematch_eve_tmin27_tmax45.mat' ;
            out = buildLookupMapTiming(mutantDir, prepend, postpend, timematfn, timematfn2, false) ;
            obj.map = out ;
        end
            
        function outstruct = findTime(obj, tfind, eps)
            %FINDTIME(tfind, eps) Find all instances with time near tfind
            %   Give the labels, folders, and time uncertainties of all
            %   stained samples matching the supplied time tfind, within
            %   a value eps of that time. Returns an output struct.
            %
            % Parameters
            % ----------
            % obj : the class instance of lookupTable
            % tfind : float or int timestamp
            % eps : optional, the allowed difference from tstamp
            
            if nargin < 3
                eps = 0.5 ;
                if nargin < 2
                    error("Must supply tfind (a time to search for) for class method findTime()")
                end
            end
            
            folders = []; 
            uncs = []; 
            dmyk = 1 ;
            % get all the labels
            labels = obj.map.keys ;
            % For each label, add struct.label, folders, uncs
            outstruct.labels = {} ;
            outstruct.folders = {} ;
            outstruct.uncs = [] ;
            for ii = 1:length(labels)
                label = labels(ii) ;
                substruct = obj.map(label{1}) ; 
                timestamps = substruct.times ;
                for jj = 1:length(timestamps)
                    tstamp = timestamps(jj) ;
                    if abs(tstamp - tfind) < eps 
                        folders{dmyk} = substruct.folders{jj} ;
                        uncs{dmyk} = substruct.uncs(jj) ;
                        dmyk = dmyk + 1 ;
                    end
                end
                % Add to the output struct
                outstruct.labels{length(outstruct.labels) + 1} = label ;
                outstruct.folders{length(outstruct.folders) + 1} = folders ;
                outstruct.uncs{length(outstruct.uncs) + 1} = uncs ;
            end
        end
        
        function outstruct = findLabel(obj, label2find)
            %FINDTIME(tfind, eps) Find all instances with given channel
            %   Give the times, folders, and time uncertainties of all
            %   stained samples matching the supplied channel 'label'
            %
            % Parameters
            % ----------
            % obj : the class instance of lookupTable
            % label : string, label name (ex 'Eve' or 'Runt')
            
            if nargin < 2
                error("Must supply label (a channel to search for) for class method findLabel()")
            end
            
            folders = []; 
            timepts = []; 
            uncs = []; 
            dmyk = 1 ;
            % get all the labels
            labels = obj.map.keys ;
            % For each label, add struct.label, folders, uncs
            outstruct.times = [] ;
            outstruct.folders = {} ;
            outstruct.uncs = [] ;
            for ii = 1:length(labels)
                label = labels(ii) ;
                if strcmp(label, label2find)
                    substruct = obj.map(label{1}) ;
                    timestamps = substruct.times ;
                    for jj = 1:length(timestamps)
                        % Add all timestamps, no matter what they are
                        timepts(dmyk) = timestamps(jj) ;
                        folders{dmyk} = substruct.folders{jj} ;
                        uncs(dmyk) = substruct.uncs(jj) ;
                        dmyk = dmyk + 1 ;
                    end
                end
            end
            
            % Add to the output struct
            outstruct.times = timepts ;
            outstruct.folders = folders ;
            outstruct.uncs = uncs ;
        end
        
        
        function outstruct = findLabelTime(obj, label2find, tfind, eps)
            %FINDTIME(tfind, eps) Find all instances with given channel
            %   Give the times, folders, and time uncertainties of all
            %   stained samples matching the supplied channel 'label'
            %
            % Parameters
            % ----------
            % obj : the class instance of lookupTable
            % label : string, label name (ex 'Eve' or 'Runt')
            % tfind : float or int timestamp
            % eps : optional, the allowed difference from tstamp
            %
            % Returns
            % -------
            % outstruct : struct
            %   The output structure with fields folders, times, uncs
            
            
            if nargin < 4
                eps = 0.5 ;
                if nargin < 3
                    error("Must supply both label and tfind for class method findLabelTime()")
                end
            end
            
            folders = []; 
            timepts = []; 
            uncs = []; 
            dmyk = 1 ;
            % get all the labels
            labels = obj.map.keys ;
            % For each label, add struct.label, folders, uncs
            outstruct.times = [] ;
            outstruct.folders = {} ;
            outstruct.uncs = [] ;
            for ii = 1:length(labels)
                label = labels(ii) ;
                if strcmp(label, label2find)
                    substruct = obj.map(label{1}) ;
                    timestamps = substruct.times ;
                    for jj = 1:length(timestamps)
                        tstamp = timestamps(jj) ;
                        if abs(tstamp - tfind) < eps 
                            timepts(dmyk) = tstamp ;
                            folders{dmyk} = substruct.folders{jj} ;
                            uncs(dmyk) = substruct.uncs(jj) ;
                            dmyk = dmyk + 1 ;
                        end
                    end
                end
            end
            
            % Add to the output struct
            outstruct.times = timepts ;
            outstruct.folders = folders ;
            outstruct.uncs = uncs ;
        end
    end
    
end

