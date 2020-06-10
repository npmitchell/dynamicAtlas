classdef dynamicAtlas < handle
    %DYNAMICATLAS Dynamic Atlas of biological images
    %   An object that handles timing datasets of dynamic 2D images 
    %
    % Example Usage
    % -------------
    % da = dynamicAtlas.dynamicAtlas('/path/to/data/', {'WT', 'TollRm9'})
    % da.buildLookup() 
    % da.makeGradientImages()
    % mapWT = da.lookup('WT') ;
    % % whose methods are
    % methods(mapWT) 
    % % To find embryos with a Runt stain
    % runts = mapWT.findLabel('Runt') 
    % % or equivalently
    % runts = mapWT.map('Runt') 
    % % To find embryos with a t=10 +/- 2 min
    % snaps = mapWT.findTime(10, 2) 
    % % To find Runt stains with a t=10 +/- 2 min
    % runtsnaps = mapWT.findLabelTime('Runt', 10, 2)
    % % Use dynamic datasets within the lookupMap to build master timeline
    % % To control how this is performed, toggle da.timeLineMethod
    % da.makeMasterTimeline('WT', 'Runt')
    % % Timestamp other data against the master timeline
    % % To control how this is performed, toggle da.timeStampMethod
    % da.timeStampStripe7('WT', 'Runt')
    %
    % To Do
    % -----
    % handle PIV-based timeline creation
    
    properties
        path
        genotypes
        lookup
        timeLineMethod
        timeStampMethod
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    methods (Abstract)
        
    end
    
    %---------------------------------------------------------------------
    % private methods
    %---------------------------------------------------------------------

    methods
        function da = dynamicAtlas(atlasPath, genotypes, Options)
            %DYNAMICATLAS Construct an instance of this class
            %   Create a dynamicAtlas instance
            %
            % Parameters
            % ----------
            % atlasPath : str
            %   the path containing each genotype directory
            % genotypes : optional cell array of strings (default = all)
            %   the genotypes to include in the atlas object
            % Options : struct with fields
            %   timeLineMethod : optional str (default = 'realspace')
            %       Method string specifier for building master timeLine(s)
            %   timeStampMethod : optional str (default = 'stripe7')
            %       Method string specifier for time stamping
            %   timerfn : str, default='timematch_*_*stripe7_chisq.mat'
            %       name of file to use to obtain timestamp for each
            %       embryo
            %   prepend : str, default='MAX_Cyl1_2_000000_c*_rot_scaled_view1' 
            %       filename search string without extension for pullback
            %   exten : str, default='.tif'
            %       filename extension for pullback
            %
            % Returns
            % -------
            % da : dynamicAtlas class instance
            %
            
            if nargin < 3
                Options = struct() ;
            end
            
            disp('Constructing dynamicAtlas')
            da.path = atlasPath ;
            
            % Add paths for methods (necessary for some MATLAB versions?)
            % da.addPaths()
            
            % Declare which genotypes to include
            if nargin > 1 && ~isempty(genotypes)
                da.genotypes = genotypes ;
            else
                % build genotype from contained directories, ignoring code
                % and extra and timing subdirs
                [~, da.genotypes] = subdirs(da.path, {'timing', 'code', 'extra'});
            end
            
            % Declare property timeLineMethod
            if nargin > 2 && isfield(Options, 'timeLineMethod')
                da.timeLineMethod = Options.timeLineMethod ;
                Options = rmfield(Options, 'timeLineMethod') ;
            else
                da.timeLineMethod = 'realspace' ;
            end
            % Declare property timeStampMethod
            if nargin > 2 && isfield(Options, 'timeStampMethod')
                da.timeStampMethod = Options.timeStampMethod ;
                Options = rmfield(Options, 'timeStampMethod') ;
            else
                da.timeStampMethod = 'stripe7' ;
            end
            
            % Populate the lookup property
            da.buildLookup(da.genotypes, Options) 
        end
        
        function addPaths(da)
            % ADDPATHS(da) Add paths needed for all routines
            %   Add paths needed for subsequent methods in dynamicAtlas
            
            % Now, instead use:
            % addpath(genpath('dynamicAtlasCode')) ;
            
            tmp = what('dynamicAtlas') ;
            tlaDir = tmp.path ;
            % splitpath = regexp(activefn, filesep, 'split') ;
            % tlaDir = '' ;
            % for qq=1:length(splitpath)
            %     tlaDir = [tlaDir splitpath{qq} filesep] ;
            % end
            
            addpath(fullfile(tlaDir, 'data_handling'))
            addpath(fullfile(tlaDir, 'tiff_handling'))
            addpath(fullfile(tlaDir, 'nanconv')) ;
            addpath(fullfile(tlaDir, 'statistics', 'chisquared')) ;
            addpath(fullfile(tlaDir, 'matchTime')) ;
            addpath(fullfile(tlaDir, 'polyparci'))
            addpath(fullfile(tlaDir, 'basics'))
            addpath(fullfile(tlaDir, 'plotting'))
            addpath(fullfile(tlaDir, 'plotting', 'ploterr'))
            addpath(fullfile(tlaDir, 'stripeExtraction'))
            
            % Add external packages
            fmDir = fullfile(tlaDir, 'external', 'toolbox_fast_marching', 'toolbox_fast_marching') ;
            addpath(fmDir)
            addpath(fullfile(fmDir, 'mex')) ;
            addpath(fullfile(fmDir, 'toolbox')) ;
            addpath(fullfile(tlaDir, 'external', 'export_fig')) ;
            addpath(fullfile(tlaDir, 'external', 'freezeColors')) ;
        end
        
        function buildLookup(da, genotypes_subset, Options)
            % BUILDLOOKUP(da)
            %   Create lookup map object for each genotype in dynamicAtlas
            % 
            % Parameters
            % ----------
            % 
            %
            da.lookup = containers.Map() ;
            if nargin < 2
                genotypes_todo = da.genotypes ;
                Options = struct() ;
            else
                genotypes_todo = genotypes_subset ;
            end
            
            
            genotypes_todo
            % Go through each genotype directory, build lookupMap class
            for kk=1:length(genotypes_todo)
                genoDir = fullfile(da.path, genotypes_todo{kk}) ;   
                % build Containers.map of lookupMaps, one for each genoDir
                lum = dynamicAtlas.lookupMap(genoDir, Options) ;
                da.lookup(da.genotypes{kk}) = lum ;
            end
        end
        
        makeGradientImages(da, sigmas, steps, cdf_minmax, overwrite)
        
        function makeMasterTimeline(da, genotype, label, Options)
            % MAKEMASTERTIMELINE(genotype, label, Options)
            %   Build a master timeline for this genotype based on the 
            %   dynamic pullbacks contained in genotype/label/
            %   
            % Parameters
            % ----------
            % da : dynamicAtlas class instance
            % genotype : str
            % stain : str
            % Options : struct with optional fields
            %   preview : optional bool, view intermediate results
            %   overwrite : optional bool, overwrite previous master timeline         
            %   thres : optional float, threshold intensity for binarization 
            %   ssfactor : optional int, subsampling factor
            %   apCijFrac : optional float between 0-1, amount of anterior end to ignore
            %
            % Returns
            % -------
            %
            % Outputs
            % -------
            % dynamicAtlas.path/timing/genotype/label/realspace_corr_ss%02d/
            % dynamicAtlas.path/timing/genotype/label/stripe7corr_ss%02d/
            % dynamicAtlas.path/timing/genotype/label/timeline_ss%02d_<corr_method>corr/
            % dynamicAtlas.path/genotype/label/embryoID/timematch_<method>_dynamic.mat
            % dynamicAtlas.path/genotype/label/embryoID/timematch_<method>_dynamic.txt

            if strcmp(da.timeLineMethod, 'realspace')
                da.makeMasterTimeLineRealspace(genotype, label, Options)
            elseif strcmp(da.timeLineMethod, 'piv')
                makeMasterTimeLinePIV(da, genotype, label, Options)
            else
                error(['dynamicAtlas.timeLineMethod not recognized: ',...
                    da.timeLineMethod])
            end
        end
        
        makeMasterTimeLineStripe7(da, genotype, label)
                
        function timeStamp(da, genotype, label, Options)
            % TIMESTAMP(da, genotype, label, Options)
            %   Find appropriate timestamp for data stained by 'label' for all fixed
            %   pullbacks compared to master timeline
            %   
            % Parameters
            % ----------
            % da : dynamicAtlas class instance
            %   the dynamicAtlas for which we find timestamps
            % genotype : str
            % stain : str
            % Options : struct with optional fields
            %   save_fancy : optional bool, default=true
            %   overwrite : optional bool, default=false
            %       overwrite previous results
            %   hands_on : optional bool, default=true
            %       use interactive domain selection for chisquare fitting
            %   cdf_min : optional float, default=0.01
            %       intensity cumulative distribution function minimum 
            %       cutoff
            %   cdf_max : optional float, default=0.999 
            %       intensity cumulative distribution function maximum 
            %       cutoff
            %   sigma : optional float, default=20 
            %       smoothing used in the stripe ID in iLastik
            %   optimize_trans : optional bool, default=true
            %       allow translation of the stripe7 curve in fitting to reference
            %       curves
            %   timeLineLabel : optional str, default=<same as label>
            %       allows the <label> fixed data to  be compared to 
            %       For ex, label=Runt (fixed data is Runt), but
            %       masterTimeLineLabel=Eve (live data is Eve).
            %   timeLineLeadingTrailing : optional str, default='leading'
            %       whether to use leading edge of timeLine stripe7 data or trailing
            %       edge
            %
            % Returns
            % -------
            %
            % Outputs
            % -------
            % dynamicAtlas.path/genotype/label/embryoID/timematch_curve7_chisq.mat
            % dynamicAtlas.path/genotype/label/embryoID/timematch_curve7_chisq.txt
            timeStampStripe7(da, genotype, label, Options)
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % QUERY METHODS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function qs = findTime(da, tfind, eps)
            %FINDTIME(tfind, eps) Find all instances with time near tfind
            %   Give the labels, folders, and time uncertainties of all
            %   stained samples matching the supplied time tfind, within
            %   a value eps of that time. Returns an output struct.
            %
            % Parameters
            % ----------
            % obj : the class instance of lookup
            % tfind : float or int timestamp
            % eps : optional, the allowed difference from tstamp
            
            if nargin < 3
                eps = 0.5 ;
                if nargin < 2
                    error("Must supply tfind (a time to search for) for class method findTime()")
                end
            end
            
            % Grab correctly timed embryos for each genotype
            genoKeys = keys(da.lookup) ;
            genotypes = {} ;
            labels = {} ;
            folders = {} ;
            names = {} ;
            embryoIDs = {} ;
            times = [] ;
            uncs = [] ;
            tiffpages = [] ;
            for ii = 1:length(genoKeys)
                genoKey = genoKeys{ii} ;
                lum = da.lookup(genoKey) ;
                tstruct = lum.findTime(tfind, eps) ;
                if ~isempty(tstruct.embryoIDs)
                    for qq = 1:length(tstruct.embryoIDs)
                        genotypes{length(genotypes) + qq} = genoKey ; 
                        labels{length(labels) + qq} = tstruct.labels{qq} ; 
                        folders{length(folders) + qq} = tstruct.folders{qq} ; 
                        names{length(names) + qq} = tstruct.names{qq} ; 
                        embryoIDs{length(embryoIDs) + qq} = tstruct.embryoIDs{qq} ; 
                    end
                    times = [times, tstruct.times] ;
                    uncs = [times, tstruct.uncs] ;
                    tiffpages = [times, tstruct.tiffpages] ;
                end 
            end
            outstruct = struct() ;
            outstruct.genotypes = genotypes ;
            outstruct.labels = labels ;
            outstruct.folders = folders ;
            outstruct.names = names ;
            outstruct.embryoIDs = embryoIDs ;
            outstruct.times = times ;
            outstruct.uncs = uncs ;
            outstruct.tiffpages = tiffpages ;
            qs = queriedSample(outstruct) ;
        end
        
        function qs = findStaticGenotypeLabel(da, genotype, label)
            %FINDDYNAMICGENOTYPELABEL(genotype, label2find) Find dynamic embryos with label
            %   Give the times, folders, and time uncertainties of all
            %   live samples matching the supplied channel 'label'
            %
            % Parameters
            % ----------
            % genotype : str, the genotype in which to search
            % label : string, label name (ex 'Eve' or 'Runt')
            % 
            % Outputs
            % -------
            % lum : struct with fields
            %   folders : Nx1 cell array of paths to label data
            %   names : Nx1 cell array of label data filenames
            %   embryoIDs : Nx1 cell array of embryoIDs
            %   times : Nx1 cell array of timestamps (each could be array)
            %   unc : Nx1 cell array of timestamp uncertainties
            %   nTimePoints : Nx1 int array of #timepoints in each dataset
            qs = da.lookup(genotype).findStaticLabel(label) ;
        end
        
        function lum = findDynamicGenotypeLabel(da, genotype, label)
            %FINDDYNAMICGENOTYPELABEL(genotype, label2find) Find dynamic embryos with label
            %   Give the times, folders, and time uncertainties of all
            %   live samples matching the supplied channel 'label'
            %
            % Parameters
            % ----------
            % genotype : str, the genotype in which to search
            % label : string, label name (ex 'Eve' or 'Runt')
            % 
            % Outputs
            % -------
            % lum : struct with fields
            %   folders : Nx1 cell array of paths to label data
            %   names : Nx1 cell array of label data filenames
            %   embryoIDs : Nx1 cell array of embryoIDs
            %   times : Nx1 cell array of timestamps (each could be array)
            %   unc : Nx1 cell array of timestamp uncertainties
            %   nTimePoints : Nx1 int array of #timepoints in each dataset
            qs = da.lookup(genotype).findDynamicLabel(label) ;
        end
            
        function qs = findGenotypeLabel(da, genotype, label)
            % FINDGENOTYPELABEL(genotype, label)
            %
            % Parameters
            % ----------
            % genotype : str, the genotype in which to search
            % label : string, label name (ex 'Eve' or 'Runt')
            %
            % Outputs
            % -------
            % lum : struct with fields
            %   folders : Nx1 cell array of paths to label data
            %   names : Nx1 cell array of label data filenames
            %   embryoIDs : Nx1 cell array of embryoIDs
            %   times : Nx1 cell array of timestamps (each could be array)
            %   unc : Nx1 cell array of timestamp uncertainties
            %   nTimePoints : Nx1 int array of #timepoints in each dataset
            %
            if strcmp(da.timeStampMethod, 'stripe7')
                qs = da.lookup(genotype).findLabel(label) ;
            else
                error(['dynamicAtlas.timeStampMethod not recognized: ',...
                    da.timeStampMethod])
            end
        end
        
        function qs = findEmbryo(da, embryoID) 
            %FINDEMBRYO(embryoID) Find all instances with given embryo
            %   Give the names, folder, times, and time uncertainties for
            %   stained samples matching the supplied embryoID.
            %   Assumes that a given embryo has only one genotype.
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
            
            % Check that embryoID is a string
            if ~isa(embryoID, 'char')
                embryoID = num2str(embryoID) ;
            end
            
            genoKeys = keys(da.lookup) ;
            found = false ;
            for ii = 1:length(genoKeys)
                if ~found
                    genoKey = genoKeys{ii} ;
                    lum = da.lookup(genoKey) ;
                    estruct = lum.findEmbryo(embryoID) ;
                    % If the embryo was part of this genotype, then declare
                    % that we found the embryo
                    if ~isempty(estruct.labels)
                        found = true ;
                    end
                end
            end
            qs = dynamicAtlas.queriedSample(estruct) ;
            
        end
    end
end

