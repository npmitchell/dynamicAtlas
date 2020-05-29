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
        function da = dynamicAtlas(atlasPath, genotypes, timeLineMethod,...
                timeStampMethod)
            %DYNAMICATLAS Construct an instance of this class
            %   Create a dynamicAtlas instance
            %
            % Parameters
            % ----------
            % atlasPath : str
            %   the path containing each genotype directory
            disp('Constructing dynamicAtlas')
            da.path = atlasPath ;
            
            % Add paths for methods
            da.addPaths()
            
            if nargin > 1
                da.genotypes = genotypes ;
            else
                % build genotype from contained directories
                [~, da.genotypes] = subdirs(da.path, {'timing', 'code', 'extra'});
            end
            % Declare property timeLineMethod
            if nargin > 2
                da.timeLineMethod = timeLineMethod ;
            else
                da.timeLineMethod = 'realspace' ;
            end
            % Declare property timeStampMethod
            if nargin > 3
                da.timeStampMethod = timeStampMethod ;
            else
                da.timeStampMethod = 'stripe7' ;
            end
        end
        
        function addPaths(da)
            % ADDPATHS(da) Add paths needed for all routines
            %   Add paths needed for subsequent methods in dynamicAtlas
            tlaDir = what('dynamicAtlas').path ;
            % splitpath = regexp(activefn, filesep, 'split') ;
            % tlaDir = '' ;
            % for qq=1:length(splitpath)
            %     tlaDir = [tlaDir splitpath{qq} filesep] ;
            % end
            addpath(fullfile(tlaDir, 'data_handling'))
            addpath(fullfile(tlaDir, 'tiff_handling'))
            addpath(fullfile(tlaDir, 'nanconv')) ;
            fmDir = fullfile(tlaDir, 'toolbox_fast_marching/toolbox_fast_marching/') ;
            addpath(fmDir)
            addpath(fullfile(fmDir, 'mex')) ;
            addpath(fullfile(fmDir, 'toolbox')) ;
            addpath(fullfile(tlaDir, 'statistics', 'chisquared')) ;
            addpath(fullfile(tlaDir, 'matchTime')) ;
            addpath(fullfile(tlaDir, 'polyparci'))
            addpath(fullfile(tlaDir, 'basics'))
            addpath(fullfile(tlaDir, 'plotting'))
            addpath(fullfile(tlaDir, 'plotting', 'ploterr'))
            addpath(fullfile(tlaDir, 'stripeExtraction'))
        end
        
        function buildLookup(da, genotypes_subset)
            % BUILDLOOKUP(da)
            %   Create lookup map object for each genotype in dynamicAtlas
            da.lookup = containers.Map() ;
            if nargin < 2
                genotypes_todo = da.genotypes ;
            else
                genotypes_todo = genotypes_subset ;
            end
            % Go through each genotype directory, build lookupMap class
            for kk=1:length(genotypes_todo)
                genoDir = fullfile(da.path, genotypes_todo{kk}) ;   
                % build Containers.map of lookupMaps, one for each genoDir
                lum = dynamicAtlas.lookupMap(genoDir) ;
                da.lookup(da.genotypes{kk}) = lum ;
            end
        end
                
        makeGradientImages(da, sigmas, steps)
        
        function makeMasterTimeline(da, genotype, label)
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
            %   preview
            %   overwrite
            %   thres : float
            %   ssfactor : int 
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
                da.makeMasterTimeLineRealspace(genotype, label)
            elseif strcmp(da.timeLineMethod, 'piv')
                makeMasterTimeLinePIV(da, genotype, label)
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
        
        function lum = findStaticGenotypeLabel(da, genotype, label)
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
            lum = da.lookup(genotype).findStaticLabel(label) ;
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
            lum = da.lookup(genotype).findDynamicLabel(label) ;
        end
            
        function lum = findGenotypeLabel(da, genotype, label)
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
                lum = da.lookup(genotype).map(label) ;
            else
                error(['dynamicAtlas.timeStampMethod not recognized: ',...
                    da.timeStampMethod])
            end
        end
        
        function estruct = findEmbryo(da, embryoID) 
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
            
        end
    end
end

