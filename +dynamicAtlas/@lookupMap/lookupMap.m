classdef lookupMap < handle
    %LOOKUPTABLE Look up timing or channels from a containers.Map object
    %   Class for a lookuptable for finding times and channels
    %
    % Example
    % -------
    % a = lookup('../WT')  % pass a genotype directory
    % a = a.buildLookup() 
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
        genoDir
        timerfn
        prepend
        exten
        map = containers.Map
        % map contains as many keys as stainings. For ex, 'Runt' is a key
        % that returns a struct with fields
        %   times : 1xN float array
        %       the timestamps for each embryo
        %   uncs : 1xN float array
        %       the uncertainty in time for each embryo
        %   folders : 1xN string cell
        %       the folder containing the pullback for each embryo
        %   names : 1xN string cell
        %       file name within folder of the pullback image
        %   embryoIDs : 1xN string cell
        %       unique date identifier for each embryo
        %   nTimePoints : 1xN int array
        %       number of timepoints in the pullback
        
    end
    
    
    methods
        function obj = lookupMap(genoDir, Options)
            % LOOKUPMAP(genoDir, Options)
            %   Instantiate a lookupMap object for genotype of genoDir
            %
            % Parameters
            % ----------
            % genoDir : str
            %   path to the genotype data directory
            % Options : optional struct with fields
            %   timerfn : str, default='timematch_curve7_chisq.mat'
            %       name of file to use to obtain timestamp for each
            %       embryo
            %   prepend : str, default='MAX_Cyl1_2_000000_c*_rot_scaled_view1' 
            %       filename search string without extension for pullback
            %   exten : str, default='.tif'
            %       filename extension for pullback
            %
            
            % name of timer filename for loading times or cell array of
            % timer filenames in rank preference order for timing
            timerfn = {'timematch_*_*stripe7_chisq.mat', ...
                'timematch_*_cephallicFurrowOnset.txt'} ;
            % filename for pullback without filetype extension
            prepend = 'MAX_Cyl1_2_00000*_c*_rot_scaled_view1' ;
            % string after prepend for the pullback data
            exten = '.tif' ;
            obj.genoDir = genoDir ;
            if nargin < 1
                error('Must supply a directory to use (genoDir) for class instantiation')
            end
            if nargin > 1
                if isfield(Options, 'timerfn')
                    timerfn = Options.timerfn ;
                end
                if isfield(Options, 'prepend')
                    prepend = Options.prepend ;
                end
                if isfield(Options, 'exten')
                    exten = Options.exten ;
                end
            end
            disp(['Options for lookupMap > ' timerfn])
            
            options = struct() ;
            if isfield(Options, 'save_map')
                options.save_map = Options.save_map ;
            end
            if isfield(Options, 'labels')
                options.labels = Options.labels ;
            end
            
            obj.timerfn = timerfn ;
            obj.prepend = prepend ;
            obj.exten = exten ;
            
            % Assign to property of lookupMap instance
            obj.map = buildLookupMap(obj, options) ;
        end
        
        %BUILDLOOKUPMAP Construct the lookup map
        out = buildLookupMap(obj, options) ;
        
        function saveLookupMap(obj)
            %SAVELOOKUPMAP Save the constructed lookup map
            map = obj.map ;
            save_map(fullfile(obj.genoDir, 'lookuptable_containersMap.mat'), 'map') ;
        end
            
        function qs = findTime(obj, tfind, eps)
            %FINDTIME(tfind, eps) Find all instances with time near tfind
            %   Give the labels, folders, and time uncertainties of all
            %   samples matching the supplied time tfind, within
            %   a value eps of that time. Returns an output struct.
            %
            % Parameters
            % ----------
            % obj : the class instance of lookup
            % tfind : float or int timestamp
            % eps : optional, the allowed difference from tstamp
            %
            % Returns
            % -------
            % qs : QueriedSample instance
            %   QueriedSample containing all samples with matching time
            
            if nargin < 3
                eps = 0.5 ;
                if nargin < 2
                    error("Must supply tfind (a time to search for) for class method findTime()")
                end
            end
            
            folders = {}; 
            embryoIDs = {} ;
            etimes = [] ;
            uncs = [] ; 
            dmyk = 1 ;
            % get all the labels
            labels = obj.map.keys ;
            % For each label, add struct.label, folders, uncs
            outstruct.labels = {} ;
            outstruct.names = {} ;
            outstruct.folders = {} ;
            outstruct.embryoIDs = {} ;
            outstruct.times = [] ;
            outstruct.uncs = [] ;
            outstruct.tiffpages = [] ;
            for ii = 1:length(labels)
                label = labels{ii} ;
                disp(['Searching ' label ' label for timestamp'])
                substruct = obj.map(label) ;
                timestamps = substruct.times ;
                for jj = 1:length(timestamps)
                    % Grab all timing + uncertainties for this embryoID
                    tstamps = timestamps{jj} ;
                    uncstamps = substruct.uncs{jj} ; % possibly multivalued
                    for kk = 1:length(tstamps)
                        % Grab the timing + uncertainty for this TIFF page
                        tstamp = tstamps(kk) ;
                        uncstamp = uncstamps(kk) ;
                        if abs(tstamp - tfind) < eps 
                            % populate the metadata for this entry
                            labels4struct{dmyk} = label ; 
                            names{dmyk} = substruct.names{jj} ; 
                            folders{dmyk} = substruct.folders{jj} ;
                            embryoIDs{dmyk} = substruct.embryoIDs{jj} ;
                            etimes(dmyk) = tstamp ;
                            uncs(dmyk) = uncstamp ;
                            tiffpages(dmyk) = kk ;
                            % update the index of the cells we are
                            % building
                            dmyk = dmyk + 1 ;
                        end
                    end
                end
            end
            
            % Add to the output struct
            if ~isempty(embryoIDs)
                outstruct.labels = labels4struct ;
                outstruct.folders = folders ;
                outstruct.names = names ;
                outstruct.embryoIDs = embryoIDs ;
                length(outstruct.times)
                outstruct.times = etimes ;
                outstruct.uncs = uncs ;
                outstruct.tiffpages = tiffpages ;
                outstruct
                qs = dynamicAtlas.queriedSample(outstruct) ;
            else
                qs = dynamicAtlas.queriedSample(struct()) ;
            end
        end
        
        function qs = findLabel(obj, label2find)
            %FINDLABEL(label2find) Find all embryos with given stain
            %   Give the times, folders, and time uncertainties of all
            %   stained samples matching the supplied channel 'label'
            %
            % Parameters
            % ----------
            % obj : the class instance of lookup
            % label : string, label name (ex 'Eve' or 'Runt')
            
            if nargin < 2
                error("Must supply label (a channel to search for) for class method findLabel()")
            end
            qs = dynamicAtlas.queriedSample(obj.map(label2find)) ;
        end
        
        function qs = findStaticLabel(obj, label2find)
            %FINDSTATICLABEL(label2find) Find dynamic embryos with label
            %   Give the times, folders, and time uncertainties of all
            %   live samples matching the supplied channel 'label'
            %
            % Parameters
            % ----------
            % obj : the class instance of lookup
            % label2find : string, label name (ex 'Eve' or 'Runt')
            qs = dynamicAtlas.queriedSample(...
                buildStructWrtTime(obj, 'static', label2find)) ;
        end
        
        function qs = findDynamicLabel(obj, label2find)
            %FINDDYNAMICLABEL(label2find) Find dynamic embryos with label
            %   Give the times, folders, and time uncertainties of all
            %   live samples matching the supplied channel 'label'
            %
            % Parameters
            % ----------
            % obj : the class instance of lookup
            % label2find : string, label name (ex 'Eve' or 'Runt')
            qs = dynamicAtlas.queriedSample(buildStructWrtTime(obj, 'dynamic', label2find)) ;
        end
        
        function qs = findStaticLabelTime(obj, label2find, tfind, deltaT)
            %FINDSTATICLABELTIME(label2find, tfind, deltaT) Find dynamic embryos with label
            %   Give the times, folders, and time uncertainties of all
            %   live samples matching the supplied channel 'label'
            %   at timestamps in range tfind-deltaT < t < tfind+deltaT.
            %
            % Parameters
            % ----------
            % obj : the class instance of lookup
            % label2find : string, label name (ex 'Eve' or 'Runt')
            % tfind : timestamp to search for (within tfind+/- deltaT)
            % deltaT : halfwidth of range of times for which to search
            %
            % Returns
            % -------
            % qs : queriedSample instance
            %   collection of the data matching criteria, as queriedSample 
            %   class instance
            qs = dynamicAtlas.queriedSample(buildStructWrtTime(obj, 'static', label2find, tfind, deltaT)) ;
        end
        
        function qs = findDynamicLabelTime(obj, label2find, tfind, deltaT)
            %FINDDYNAMICLABELTIME(label2find, tfind, deltaT) Find dynamic embryos with label
            %   Give the times, folders, and time uncertainties of all
            %   dynamic/live samples matching the supplied channel 'label' 
            %   at timestamps in range tfind-deltaT < t < tfind+deltaT.
            %
            % Parameters
            % ----------
            % obj : the class instance of lookup
            % label2find : string, label name (ex 'Eve' or 'Runt')
            % tfind : timestamp to search for (within tfind+/- deltaT)
            % deltaT : halfwidth of range of times for which to search
            %
            % Returns
            % -------
            % qs : queriedSample instance
            %   collection of the data matching criteria, as queriedSample 
            %   class instance
            if nargin < 3
                tfind = 1 ;
            end
            if nargin < 4
                deltaT = 1 ;
            end
            tmp = buildStructWrtTime(obj, 'dynamic', label2find, tfind, deltaT) ;
            
            qs = dynamicAtlas.queriedSample(tmp) ;
        end
               
        function qs = findLabelTime(obj, label2find, tfind, eps)
            %FINDLABELTIME(label2find, tfind, eps)
            %   Give the times, folders, and time uncertainties of all
            %   stained samples matching the supplied channel 'label'
            %
            % Parameters
            % ----------
            % obj : the class instance of lookup
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
            
            folders = {}; 
            names = {} ;
            labels4struct = {} ;
            embryoIDs = {}; 
            timepts = []; 
            uncs = []; 
            tiffpages = [] ;
            nTimePoints = [] ;
            dmyk = 1 ;
            % get all the labels
            labels = obj.map.keys ;
            % For each label, add struct.label, folders, uncs
            for ii = 1:length(labels)
                label = labels(ii) ;
                if strcmp(label, label2find)
                    substruct = obj.map(label{1}) ;
                    timestamps = substruct.times ;
                    % Consider each embryo and look for time matches
                    for jj = 1:length(timestamps)
                        tstamps = timestamps{jj} ;
                        uncstamps = substruct.uncs{jj} ;
                        % go through each page of the tiff and look for
                        % time matches
                        for kk = 1:length(tstamps) 
                            tstamp = tstamps(kk) ;
                            uncstamp = uncstamps(kk) ;
                            if abs(tstamp - tfind) < eps
                                labels4struct{dmyk} = label2find ;
                                folders{dmyk} = substruct.folders{jj} ;
                                names{dmyk} = substruct.names{jj} ;
                                embryoIDs{dmyk} = substruct.embryoIDs{jj} ;
                                timepts(dmyk) = tstamp ;
                                uncs(dmyk) = uncstamp ;
                                tiffpages(dmyk) = kk ;
                                nTimePoints(dmyk) = substruct.nTimePoints(jj) ;
                                dmyk = dmyk + 1 ;
                            end
                        end
                    end
                end
            end
            
            % Add to the output struct
            outstruct.labels = labels4struct ;
            outstruct.folders = folders ;
            outstruct.names = names ;
            outstruct.embryoIDs = embryoIDs ;
            outstruct.times = timepts ;
            outstruct.uncs = uncs ;
            outstruct.tiffpages = tiffpages ;
            outstruct.nTimePoints = nTimePoints ;
            
            qs = dynamicAtlas.queriedSample(outstruct) ;
        end
        
        function estruct = findEmbryo(obj, embryoID) 
            %FINDEMBRYO(embryoID) Find all instances with given embryo
            %   Give the names, folder, times, and time uncertainties for
            %   stained samples matching the supplied embryoID
            %
            % Parameters
            % ----------
            % obj : the class instance of lookup
            % embryoID : string, embryo identification name 
            %           (ex '201904011200')
            %
            % Returns
            % -------
            % estruct : struct
            %   The output structure with fields labels, folders, times, uncs
            
            elabels = {} ;
            folders = {} ; 
            names = {} ;
            embryoIDs = {} ;
            timepts = {} ;  % contains arrays for dynamic data, float for fixed data
            uncs = {} ; 
            nTimePoints = [] ;
            dmyk = 1 ;
            % get all the labels
            labels = obj.map.keys ;
            % For each label, add struct.label, folders, uncs
            for ii = 1:length(labels)
                label = labels(ii) ;
                label = label{1} ;
                substruct = obj.map(label) ;
                eIDs = substruct.embryoIDs ;
                for jj = 1:length(eIDs)
                    if strcmp(eIDs{jj}, embryoID)
                        elabels{dmyk} = label ;
                        folders{dmyk} = substruct.folders{jj} ;
                        names{dmyk} = substruct.names{jj} ;
                        embryoIDs{dmyk} = substruct.embryoIDs{jj} ;
                        timepts{dmyk} = substruct.times{jj} ;
                        uncs{dmyk} = substruct.uncs{jj} ;
                        nTimePoints(dmyk) = substruct.nTimePoints(jj) ;
                        dmyk = dmyk + 1 ;
                    end
                end
            end
            
            % Add to the output struct
            estruct.labels = elabels ;
            estruct.folders = folders ;
            estruct.names = names ;
            estruct.embryoIDs = embryoIDs ;
            estruct.times = timepts ;
            estruct.uncs = uncs ;
            estruct.nTimePoints = nTimePoints ; 
            qs = dynamicAtlas.queriedSample(estruct) ;
        end
        
        % Util methods 
        outstruct = buildStructWrtTime(obj, timestr, label2find, tfind, deltaT) ;
        
    end
    
end

