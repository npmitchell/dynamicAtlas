classdef queriedSample < handle
    %QUERIEDSAMPLE Object for manipulating subset of atlas
    %   Class for manipulating a queried subset of the dynamicAtlas
    %   Holds a subset of embryos in the atlas and offers methods for
    %   interrogating the properties of those embryos, such as loading the
    %   data, smoothed data, gradients, flow fields.
    %
    %   Note that a queriedSample can slice across embryos, labels, 
    %   genotypes, or timestamps. There are few restrictions about what
    %   kinds of embryos can be contained in the queried sample. For
    %   example, da.findTime() will return a queriedSample of all embryos
    %   of any genotype in the da (dynamicAtlas instance) that match a
    %   timestamp range.
    %
    %   Indexing convention is always Ntps x spaceDim1 x spaceDim2.
    %   
    % Properties
    % ----------
    % meta: struct with fields
    %     folders : Nx1 cell array of paths to label data
    %     names : Nx1 cell array of label data filenames
    %     embryoIDs : Nx1 cell array of embryoIDs
    %     times : Nx1 cell array of timestamps (each could be array)
    %     unc : Nx1 cell array of timestamp uncertainties
    %     nTimePoints : Nx1 int array of #timepoints in each dataset
    % data : initially empty cell of data images
    %     can be populated via qs.getData()
    % gradients : initially empty cell of gradient images
    %     can be populated via qs.getGradients()
    % smooth : initially empty cell of smoothed data images
    %     can be populated via qs.getSmooth()
    % piv : initially empty cell of velocity fields
    %     can be populated via qs.getPIV()
    %
    % Example
    % -------
    % da = dynamicAtlas('./', {'WT'}) 
    % qs = da.findGenotypeLabel('WT', 'Eve')
    % qs.getData()
    % qs.getGradients('dx', 10)
    % qs.getSmooth()
    % qs.buildPIVStack(da, genotype, labels, deltaTime)
    %
    % NPMitchell 2020
    
    properties
        meta                    % struct with fields labels, folders, names, ...
                                %       embryoIDs, times, uncs, tiffpages, nTimePoints
        data                    % cell array of images
        gradients = struct()    % struct with fields dx, dy, mag
        smooth                  % cell array of images
        dynamicAtlas            % dynamicAtlas object from which qs is made
        piv = struct('vx', [], ...
            'vy', [], ...
            'rescaleFactor', 0.4)   % cell array of velocity fields
        meanPIV = struct('vx', [], ...
            'vy', [])           % cell array of average velocity field across sample
    end
    
    
    methods
        function obj = queriedSample(sampleStruct)
            % sampleStruct:  struct with fields 
            %   labels, folders, names, ...
            %   embryoIDs, times, uncs, tiffpages, nTimePoints
            obj.meta = sampleStruct ;
            
            % if we want to get the common string name of all tiffs:
            % if all(strcmp(Scell(:), Scell{1}))
            %     obj.prependFileName = obj.meta.names{1}(all(~diff(char(Scell(:))))) ;
            % else
            %     obj.prependFileName = obj.meta.names{1}(all(~diff(char(Scell(:))))) ;
            % end
        end
        
        % function obj = findTime()
        % It is not idiomatic to findTime() here, instead use LUM
        
        timeStamp(obj, da, genotype, label, options)
        
        function dataCell = getData(obj)
            %GETDATA(obj) Load the TIFFs for all data in this queriedSample
            %
            % Returns
            % -------
            % dataCell : N x 1 cell of 2d uint16 arrays
            %   The requested raw data images
            %
            if isempty(obj.data)
                dataCell = cell(1, length(obj.meta.folders)) ;
                for qq = 1:length(obj.meta.folders)
                    tiff_fn = fullfile(obj.meta.folders{qq}, obj.meta.names{qq}) ;


                    if isfield(obj.meta, 'tiffpages') && ~iscell(obj.meta.tiffpages)

                        %MODIFIED 2024/12/11, the if statement did nothing, put
                        %inside this if statement
                        thesePages = obj.meta.tiffpages{qq} ;

                        disp(['Loading TIFF page ', ...
                            num2str(thesePages), ...
                            ' from: ', tiff_fn])
                        try
                            dataCell{qq} = imread(tiff_fn, thesePages) ;
                        catch
                            disp(['Could not load file. Incorrectly named? '  tiff_fn])
                            dataCell{qq} = imread(tiff_fn, thesePages) ;
                        end
                    else
                        disp(['Loading TIFF stack: ' tiff_fn])
                        dataCell{qq} = loadtiff(tiff_fn) ;
                    end
                end
                obj.data = dataCell ;
            else
                % If we have requested the data to be returned as a cell,
                % grab that cell 
                if nargout > 0
                    dataCell = obj.data ;
                end
            end
        end
        
        function dataCell = getGradients(obj, xymag, sigma, step)
            %GETGRADIENTS(obj) 
            %   Load the gradient TIFFs for all data in this queriedSample
            %
            % Parameters
            % ----------
            % xymag : optional str ('x', 'y', 'mag'), default='mag'
            %   whether to grab AP gradient (x), DV gradient (y), or
            %   normed magnitude of gradient
            % sigma : optional int, default=10
            %   length scale over which data is smoothed before gradient
            %   computed
            % step : optional int, default=1
            %   length scale over which gradient is computed
            %
            % Returns
            % -------
            % dataCell : N x 1 cell of 2d uint16 arrays
            %   The requested gradient images
            %
            
            % Default args
            if nargin < 4
                step = 1;
            end
            if nargin < 3
                sigma = 10 ;
            end
            if nargin < 2
                xymag = 'mag' ;
            end
            
            % Grab data
            if isempty(obj.data)
                dataCell = cell(1, length(obj.meta.folders)) ;
                for qq = 1:length(obj.meta.folders)
                    gradsubdir = sprintf('sigma%03d_step%03d', sigma, step) ;
                    if contains(xymag, 'x')
                        insertion = 'dx_' ;
                    elseif contains(xymag, 'y')
                        insertion = 'dy_' ;
                    else
                        insertion = '' ;
                    end
                    subsubdir = ['gradient_' insertion 'magnitude'] ;
                    name = [obj.meta.labels{qq}, '_gradient_', ...
                        insertion, 'magnitude_', gradsubdir, '.tif'] ;
                    tiff_fn = fullfile(obj.meta.folders{qq}, ...
                        gradsubdir, subsubdir, name) ;
                    
                    % Now load the gradient image 
                    if isfield(obj.meta, 'tiffpages')
                        disp(['Loading TIFF page ', ...
                            num2str(obj.meta.tiffpages(qq)), ...
                            ' from: ', tiff_fn])
                        dataCell{qq} = imread(tiff_fn, obj.meta.tiffpages(qq)) ;
                    else
                        disp(['Loading TIFF stack: ' tiff_fn])
                        dataCell{qq} = loadtiff(tiff_fn) ;
                    end
                end
                if contains(xymag, 'x')
                    obj.gradients.dx = dataCell ;
                elseif contains(xymag, 'y')
                    obj.gradients.dy = dataCell ;
                else
                    obj.gradients.mag = dataCell ;
                end
            else
                % If we have requested the data to be returned as a cell,
                % grab that cell 
                if nargout > 0
                    dataCell = obj.data ;
                end
            end
        end
        
        function meanIm = getMeanData(obj, normalizeEach, preview)
            % Load all data and take the mean of all data in the sampling 
            %
            % Parameters
            % ----------
            % obj : queriedSample instance (self)
            % normalizeEach : bool (default = true)
            %   normalize the intensity of each frame as we sum them
            %   (recommended)
            % preview : bool (default = false)
            %   preview the averaging of all members of the queriedSample,
            %   element by element as we sum them
            %
            % Returns
            % -------
            % meanIm : 2D double array
            %   mean image intensity for all samples in queriedSample
            %   instance
            %
            if nargin < 2
                normalizeEach = true ;
            end
            if nargin < 3
                preview = false ;
            end
            dataCell = obj.getData() ;
            if isempty(dataCell)
                disp('WARNING: No data in qs! Returning empty array.')
                meanIm = [] ;
            else
                % First image
                if normalizeEach
                    toAdd = dataCell{1} ;
                    meanIm = double(toAdd) / double(max(toAdd(:))) ;
                else
                    meanIm = double(dataCell{1}) ;
                end
                % check the first image
                if preview
                    imagesc(meanIm); 
                    colorbar ;
                    pause(0.1)
                end
                
                % Add all other images and normalize at the end
                for qq = 2:length(dataCell)
                    
                    % Normalize the image or just sum it up
                    if normalizeEach
                        toAdd = dataCell{qq} ;
                        meanIm = meanIm + double(toAdd) / double(max(toAdd(:))) ;
                    else
                        meanIm = meanIm + double(dataCell{qq}) ;
                    end
                    
                    % check the current sum (normalized)
                    if preview
                        imagesc(meanIm / qq);
                        pause(0.1)
                    end
                end            
                
                % Normalize by the number of frames summed
                meanIm = double(meanIm) ./ length(dataCell) ;
            end
        end
                        
        function ensurePIV(obj, options)
            % Make sure each embryo in the queriedSample has PIV and
            % PIV_filtered measurements: one .mat per timepoint
            nEmbryos = length(obj.meta.folders) ;
            method = 'default' ;
            rescaleFactor = obj.piv.rescaleFactor ;
            if nargin < 2
                options = struct() ;
            end
            % default options
            if isfield(options, 'isf')
                if options.isf ~= rescaleFactor 
                    error('rescaling factor (isf) passed to ensurePIV does not match dataset global rescaleFactor')
                end
            else
                options.isf = rescaleFactor ;
            end
            if isfield(options, 'method')
                method = options.method ;
            end
            
            if ~isfield(options, 'overwrite')
                options.overwrite = false ;
            end
            if ~isfield(options, 'method')
                options.method = method ;
            end
            if ~isfield(options, 'overwrite')
                options.overwrite = false ;
            end
            if isfield(options, 'median_order')
                median_order = options.median_order ;
            else
                median_order = 3 ;
            end
            for qq = 1:nEmbryos

                %MODIFIED BY VJS 2024/12/29
                %MODIFIED TO SKIP DATASETS WHICH ARE FIXED
                if (obj.meta.nTimePoints(qq) == 1)
                    continue
                end
                
                % check if PIV has already been computed or if we are
                % overwriting the existing results
                if strcmpi(method, 'pivlab')
                    input_path = fullfile(obj.meta.folders{qq}, 'PIVlab') ;
                else
                    input_path = fullfile(obj.meta.folders{qq}, 'PIV') ;
                end
                fns = dir(fullfile(input_path, 'Velo*.mat')) ;
                recompute = options.overwrite || ~exist(input_path, 'dir') || ...
                    (length(fns)<(obj.meta.nTimePoints(qq)-1)) ;
                if recompute
                    options.Name = obj.meta.names{qq} ;
                    PIVTimeseries(obj.meta.folders{qq}, options) 
                end
                
                % check if filtered PIV already exists
                if strcmpi(method, 'default')
                    output_path = fullfile(obj.meta.folders{qq}, 'PIV_filtered') ;
                elseif strcmpi(method, 'pivlab')
                    output_path = fullfile(obj.meta.folders{qq}, 'PIVlab_filtered') ;
                end
                ffns = dir(fullfile(output_path, 'Velo*.mat')) ;
                if options.overwrite || ~exist(output_path, 'dir') || length(ffns) < length(fns)
                    disp('Running median filter on PIV...')
                    runMedianFilterOnPIVField(input_path, output_path, median_order);
                end
            end

            disp('Done with ensurePIV')
        end
        
        function dataCell = getSmooth(obj, sigma, step)
            %GETGRADIENTS(obj) 
            %   Load the smoothed TIFFs for all data in this queriedSample
            %
            % Parameters
            % ----------
            % sigma : optional int, default=10
            %   length scale over which data is smoothed before gradient
            %   computed
            % step : optional int, default=1
            %   length scale over which gradient is computed
            %
            % Returns
            % -------
            % dataCell : N x 1 cell of 2d uint16 arrays
            %   The requested gradient images
            %
            
            % Default args
            if nargin < 3
                step = 1;
            end
            if nargin < 2
                sigma = 10 ;
            end
            
            % Grab data
            if isempty(obj.data)
                dataCell = cell(1, length(obj.meta.folders)) ;
                for qq = 1:length(obj.meta.folders)
                    ss_subdir = sprintf('sigma%03d_step%03d', sigma, step) ;
                    name = [obj.meta.labels{qq}, '_smooth_', ...
                        ss_subdir, '.tif'] ;
                    tiff_fn = fullfile(obj.meta.folders{qq}, ...
                        ss_subdir, 'smooth', name) ;
                    
                    % Now load the smoothed image 
                    if isfield(obj.meta, 'tiffpages')
                        disp(['Loading TIFF page ', ...
                            num2str(obj.meta.tiffpages(qq)), ...
                            ' from: ', tiff_fn])
                        dataCell{qq} = imread(tiff_fn, obj.meta.tiffpages(qq)) ;
                    else
                        disp(['Loading TIFF stack: ' tiff_fn])
                        dataCell{qq} = loadtiff(tiff_fn) ;
                    end
                end
                obj.smooth = dataCell ;
            else
                % If we have requested the data to be returned as a cell,
                % grab that cell 
                if nargout > 0
                    dataCell = obj.data ;
                end
            end
        end
        
        function obj = removeEntries(obj, indicesToRemove)
            % removeEntries(queriedSample, indicesToRemove)
            % Remove some entries of the sample from the queriedSample
            % instance.
            %
            % Example Usage 
            % -------------
            % qs = qs.removeEntries(indicesToRemove)
            %
            nEntries = length(obj.meta.folders) ;
            assert(length(obj.meta.times) == nEntries)
            assert(length(obj.meta.uncs) == nEntries)
            assert(length(obj.meta.folders) == nEntries)
            assert(length(obj.meta.nTimePoints) == nEntries)
            assert(length(obj.meta.names) == nEntries)
            assert(length(obj.meta.embryoIDs) == nEntries)
            allIndices = 1:nEntries ;
            keep = setdiff(allIndices, indicesToRemove);
            obj.meta.folders = obj.meta.folders(keep) ;
            obj.meta.names = obj.meta.names(keep) ;
            obj.meta.times = obj.meta.times(keep) ;
            obj.meta.uncs = obj.meta.uncs(keep) ;
            obj.meta.nTimePoints = obj.meta.nTimePoints(keep) ;
            obj.meta.embryoIDs = obj.meta.embryoIDs(keep) ;
            
            try 
                assert(length(fields(obj.meta)) == 6)
            catch
                error('We have added new fields to queriedSample since writing this function. Trim those here')
            end
            
            % reset other fields -- todo: remove only the entries from the
            % other fields
            obj.data = [] ;
            obj.gradients = struct() ;
            obj.smooth = [] ;
            obj.piv = [] ;
        end
                
        function minTime = getMinTime(obj)
            % Measure the samllest timestamp value in the queried sample
            nEmbryos = length(obj.meta.folders) ;
            if isa(obj.meta.times, 'cell') 
                minTime = Inf ;
                for qq = 1:nEmbryos
                    minTime = min(min(obj.meta.times{qq}), minTime) ;
                end
            else
                minTime = min(obj.meta.times) ;
            end
        end
        
        function maxTime = getMaxTime(obj)
            % Measure the largest timestamp value in the queried sample
            nEmbryos = length(obj.meta.folders) ;
            if isa(obj.meta.times, 'cell') 
                maxTime = -Inf ;
                for qq = 1:nEmbryos
                    maxTime = max(max(obj.meta.times{qq}), maxTime) ;
                end
            else
                maxTime = min(obj.meta.times) ;
            end
        end
        
        function piv = getPIV(obj, options)
            % piv = getPIV(obj, options)
            % 
            % Parameters
            % ----------
            % options : optional struct with optional fields
            %   method : string specifier ('default' or 'pivlab')
            %       how the piv was computed.
            % 
            % Returns
            % -------
            % piv : struct with fields
            %   vx : cell of NtpsxPxQ double arrays
            %       velocity in x direction in pullbacks
            %   vy : cell of NtpsxPxQ double arrays
            %       velocity in y direction in pullbacks
            %   rescaleFactor : float
            %   X0 : QxP double
            %       piv evaluation x coordinates
            %   Y0 : QxP double
            %       piv evaluation y coordinates
            method = 'pivlab';
            if nargin < 2 
                options = struct() ;
            end
            if isfield(options, 'method')
                method = options.method ;
            end
            % Load PIV results for all entries in queriedSample collection
            nEmbryos = length(obj.meta.folders) ;
            rescaleFactor = obj.piv.rescaleFactor ;
            for qq = 1:nEmbryos
                if isa(obj.meta.times, 'cell')
                    % Load each in cell of timestamps
                    ntimestamps = obj.meta.nTimePoints(qq) ;
                    timestamps = 1:ntimestamps ;
                    if strcmpi(method, 'default')
                        vxcollection = zeros(46, 54, length(timestamps)) ;
                        vycollection = zeros(46, 54, length(timestamps)) ;
                        
                        % get evaluation coordinates
                        % pivCoordSys: 46 x 54
                        EdgeLength = 15;
                        % resized dimensions of piv grid --> nearly (szX_orig, szY_orig) * isf
                        rescaleFactor = 0.4 ;
                        szX = 696 ;
                        szY = 820 ;
                        % PIV evaluation coordinates in resized pixels
                        [X0,Y0] = meshgrid(EdgeLength/2:EdgeLength:(szX-EdgeLength/2), ...
                            EdgeLength/2:EdgeLength:(szY-EdgeLength/2)); 

                        obj.piv.X0 = X0 ;
                        obj.piv.Y0 = Y0 ;
                    else
                        disp('preparing for loading PIVlab results...')
                        vxcollection = zeros(length(timestamps), 86, 101) ;
                        vycollection = zeros(length(timestamps), 86, 101) ;
                        first = true ;
                    end
                    
                    for pp = 1:(ntimestamps-1)
                        try
                            % load dt, timestep between frames in minutes
                            dt = dlmread(fullfile(obj.meta.folders{qq}, 'dt.txt')) ;
                        catch
                            % default dt is 1 minute
                            dt = 1;
                        end
                        tiffpage = timestamps(pp) ;
                        try
                            if strcmpi(method, 'default')
                                pivfn = fullfile(obj.meta.folders{qq}, 'PIV_filtered', ...
                                    sprintf('VeloT_medfilt_%06d.mat', tiffpage)) ;
                                tmp = load(pivfn) ;
                                vxcollection(pp, :, :) = tmp.VX / (rescaleFactor *dt) ;
                                vycollection(pp, :, :) = tmp.VY / (rescaleFactor *dt) ;
                            else
                                pivfnRaw = fullfile(obj.meta.folders{qq}, 'PIVlab', ...
                                    sprintf('VeloT_fine_%06d.mat', tiffpage)) ;
                                pivfn = fullfile(obj.meta.folders{qq}, 'PIVlab_filtered', ...
                                    sprintf('VeloT_medfilt_fine_%06d.mat', tiffpage)) ;
                                tmp = load(pivfn) ;
                                if first
                                    tmp = load(pivfnRaw, 'convert_to_pix_per_min', 'X0', 'Y0') ;
                                    obj.piv.X0 = tmp.X0 ;
                                    obj.piv.Y0 = tmp.Y0 ;
                                    try
                                        conversion = tmp.convert_to_pix_per_min ;
                                    catch
                                        tmp = load(pivfnRaw, 'convert_to_um_per_min') ;
                                        conversion = tmp.convert_to_um_per_min / 0.2619 ;  % um/min * original size pix / um
                                    end
                                    first = false ;
                                else
                                    try
                                        tmp = load(pivfnRaw, 'convert_to_pix_per_min') ;
                                        conversion = tmp.convert_to_pix_per_min ;
                                    catch 
                                        tmp = load(pivfnRaw, 'convert_to_pix_per_min') ;
                                        conversion = tmp.convert_to_um_per_min / 0.2619;  % um/min * original size pix / um
                                    end
                                end
                                vxcollection(pp, :, :) = tmp.VX * conversion ;
                                vycollection(pp, :, :) = tmp.VY * conversion ;
                                obj.piv.units = 'pixels per min' ;
                            end
                        catch
                            disp('Warning: some timepoints have no velocity, returned these as NaN')

                            vxcollection(pp, :, :) = NaN ;
                            vycollection(pp, :, :) = NaN ;
                        end
                    end
                    obj.piv.vx{qq} = vxcollection ;
                    obj.piv.vy{qq} = vycollection ;
                else
                    % There is only one timepoint in this queriedSample?
                    timestamp = obj.meta.tiffpages(qq) ;
                    
                    if strcmpi(method, 'default')
                        pivfn = fullfile(obj.meta.folders{qq}, 'PIV_filtered', ...
                            sprintf('VeloT_medfilt_%06d.mat', timestamp)) ;
                        tmp = load(pivfn) ;
                        obj.piv.vx{qq} = tmp.VX / rescaleFactor ;
                        obj.piv.vy{qq} = tmp.VY / rescaleFactor ;
                    else
                        error('handle case here')
                    end
                end
            end
            if nargout > 0
                piv = obj.piv ;
            end
        end
        
        function meanPIV = getMeanPIV(obj, options)
            %meanPIV = getMeanPIV(obj, options)
            % Load PIV results for all entries in queriedSample collection
            %
            % Parameters
            % ----------
            % options : struct with fields
            %   method : 'default' or 'pivlab'
            % 
            % Returns
            % -------
            % meanPIV : struct with fields
            %   mean flow across all samples in the queriedSample, in units
            %   of pixels per minute. 
            %   struct with fields
            %       vx : cell of PxQxNtps double arrays
            %           velocity in x direction in pullbacks
            %       vy : cell of PxQxNtps double arrays
            %           velocity in y direction in pullbacks
            %       X0 : QxP double
            %           piv evaluation x coordinates
            %       Y0 : QxP double
            %           piv evaluation y coordinates
            
            % default options
            method = 'pivlab_filtered' ;
            preview = false ;
            if nargin < 2
                options = struct() ;
            end
            if isfield(options, 'method')
                method = options.method ;
            end
            if isfield(options, 'preview')
                preview = options.preview ;
            end
            
            nEmbryos = length(obj.meta.folders) ;

            % figure out the total number of frames in the queriedSample 
            nPages = 0 ;
            for qq = 1:nEmbryos
                try
                    nPages = nPages + length(obj.meta.tiffpages{qq}) ;
                catch
                    nPages = nPages + length(obj.meta.tiffpages(qq)) ;
                end
            end
            
            % preallocate the velocity information
            if strcmpi(method, 'default')
                vx = zeros(nPages, 46, 54) ;
                vy = zeros(nPages, 46, 54) ;
            else
                vx = zeros(nPages, 86, 101) ;
                vy = zeros(nPages, 86, 101) ;
            end
            
            rescaleFactor = obj.piv.rescaleFactor ;
            dmyk = 1 ;
            for qq = 1:nEmbryos
                dt = dlmread(fullfile(obj.meta.folders{qq}, 'dt.txt')) ;
                % if isa(obj.meta.tiffpages, 'cell') && length(obj.meta.tiffpages) > 1
                %     error('handle this case of collections here -- see above in getPIV()')
                for tpageID = 1:length(obj.meta.tiffpages{qq})
                    timestamp = obj.meta.tiffpages{qq}(tpageID); 
                    if strcmpi(method, 'default')
                        pivfn = fullfile(obj.meta.folders{qq}, 'PIV_filtered', ...
                            sprintf('VeloT_medfilt_%06d.mat', timestamp)) ;
                        if exist(pivfn, 'file')
                            tmp = load(pivfn) ;
                            vx(dmyk, :, :) = tmp.VX / (rescaleFactor * dt) ;
                            vy(dmyk, :, :) = tmp.VY / (rescaleFactor * dt) ;
                        else
                            if timestamp == obj.meta.nTimePoints(qq)
                                disp('Skipping final timepoint from meanPIV')
                                vx(dmyk, :, :) = NaN ;
                                vy(dmyk, :, :) = NaN ;
                            else
                                disp(['PIV does not exist but tp is not final tp: ' pivfn])
                                vx(dmyk, :, :) = NaN ; 
                                vy(dmyk, :, :) = NaN ;
                            end
                        end
                    elseif strcmpi(method, 'pivlab') || strcmpi(method, 'pivlab_filtered') 
                        if strcmpi(method, 'pivlab_filtered') 
                            pivfnRaw = fullfile(obj.meta.folders{qq}, 'PIVlab', ...
                                sprintf('VeloT_fine_%06d.mat', timestamp)) ;
                            pivfn = fullfile(obj.meta.folders{qq}, 'PIVlab_filtered', ...
                                sprintf('VeloT_medfilt_fine_%06d.mat', timestamp)) ;
                        else
                            pivfn = fullfile(obj.meta.folders{qq}, 'PIVlab', ...
                                sprintf('VeloT_fine_%06d.mat', timestamp)) ;
                        end
                        if exist(pivfn, 'file')
                            
                            % load conversion factor to pix per min
                            if strcmpi(method, 'pivlab_filtered') 
                                tmp = load(pivfn) ;
                                try
                                    rawPIV = load(pivfnRaw, 'convert_to_pix_per_min') ;
                                    tmp.convert_to_pix_per_min = rawPIV.convert_to_pix_per_min ;
                                catch
                                    rawPIV = load(pivfnRaw, 'convert_to_um_per_min') ;
                                    tmp.convert_to_pix_per_min = rawPIV.convert_to_um_per_min / 0.2619 ; 
                                end
                            else
                                tmp = load(pivfn) ;
                            end
                            
                            tmp = load(pivfn) ;
                            vx(dmyk, :, :) = tmp.VX / tmp.convert_to_pix_per_min ; % pix / min
                            vy(dmyk, :, :) = tmp.VY / tmp.convert_to_pix_per_min  ; % pix / min
                            assert(any(~isnan(vx(:))))
                        else
                            if timestamp == obj.meta.nTimePoints(qq)
                                disp('Skipping final timepoint from meanPIV')
                                vx(dmyk, :, :) = NaN ;
                                vy(dmyk, :, :) = NaN ;
                            else
                                disp(['PIV does not exist but tp is not final tp: ' pivfn])
                                vx(dmyk, :, :) = NaN ;
                                vy(dmyk, :, :) = NaN ;
                            end
                        end
                    end
                    dmyk = dmyk + 1 ;
                end
            end
            
            % Get evaluation points
            if strcmpi(method, 'pivlab')
                X0 = tmp.X0 ;
                Y0 = tmp.Y0 ;
                obj.meanPIV.units = 'pix per minute' ;
            else
                EdgeLength = 15;
                % resized dimensions of piv grid --> nearly (szX_orig, szY_orig) * isf
                if rescaleFactor == 0.4
                    szX = 696 ;
                    szY = 820 ;
                    % PIV evaluation coordinates in resized pixels
                    [X0,Y0] = meshgrid(EdgeLength/2:EdgeLength:(szX-EdgeLength/2), ...
                        EdgeLength/2:EdgeLength:(szY-EdgeLength/2)); 
                    X0 = X0' ;
                    Y0 = Y0' ;
                else
                    error('handle this rescaleFactor here')
                end
                obj.meanPIV.units = 'original pix per minute' ;
            end
            
            % Take mean
            obj.meanPIV.vx = squeeze(mean(vx, 1, 'omitnan')) ;
            obj.meanPIV.vy = squeeze(mean(vy, 1, 'omitnan')) ;
            obj.meanPIV.X0 = X0 ;
            obj.meanPIV.Y0 = Y0 ;
            if nargout > 0
                meanPIV = obj.meanPIV ;
            end
                       
            % Preview flow timecourse & output
            if preview
                subplot(1, 2, 2)
                quiver(X0, Y0, obj.meanPIV.vx * 10, obj.meanPIV.vy * 10, 0)
                axis equal 
                axis tight
                
                for dmyk = 1:nPages
                    subplot(1, 2, 1)
                    quiver(X0, Y0, squeeze(vx(dmyk, :, :)) * 10, ...
                        squeeze(vy(dmyk, :, :)) * 10, 0)
                    title(['page = ' num2str(dmyk)])
                    axis equal 
                    axis tight
                    pause(0.01)
                end
            end
            
        end
        
        function pivStack = buildPIVStack(obj, da, genotype, labels, deltaTime, method)
            % Build a stack of PIV flow fields over time
            %
            % Parameters
            % ----------
            % da : dynamicAtlas instance queried by current queriedSample 
            %      instance
            % genotype : genotype of the queried sample to query for
            %            velocities
            % labels : 
            % deltaTime
            % options : struct with fields
            %   method : 'default' or 'pivlab'
            %
            % Returns
            % -------
            % pivStack
            %
            
            if nargin < 5
                % half-width of search window for each timepoint
                deltaTime = 0.5 ;
            end
            if nargin < 6
                options = struct() ;
            end
            
            method = 'default' ;
            if isfield(options, 'method')
                method = options.method ;
            end
            
            if strcmpi(method, 'default')
                % Built-in pullback image size, hard-coded here:
                % ----------------------------------------------
                % Original image is 2000x2000, then resized is 0.4*(Lx,Ly), then PIV is 15x
                % smaller than the resized image. 
                % fullsizeCoordSys: 1738x2050
                % resizeCoordSys: 696 x 820
                % pivCoordSys: 46 x 54
                EdgeLength = 15;
                % rescaleFactor = 0.4;
                % szX_orig = 1738 ;
                % szY_orig = 2050 ;

                % resized dimensions of piv grid --> nearly (szX_orig, szY_orig) * isf
                if obj.piv.rescaleFactor == 0.4
                    szX = 696 ;
                    szY = 820 ;
                else
                    error('code for this rescaleFactor convention here')
                end

                % in resized pixels
                [X0,Y0] = meshgrid(EdgeLength/2:EdgeLength:(szX-EdgeLength/2), ...
                    EdgeLength/2:EdgeLength:(szY-EdgeLength/2)); 

                mint = obj.getMinTime() ;
                maxt = obj.getMaxTime() ;
                ntps = round(maxt-mint) ;
                pivStack = struct('vx', zeros(ntps, size(X0, 2), size(X0, 1)), ...
                    'vy', zeros(ntps, size(X0, 2), size(X0, 1))) ;
            else
                error('code for pivlab here')
            end
            
            tidx = 1 ;
            for tt = mint:maxt
                disp(['obtaining mean PIV for t=' num2str(tt)])
                snap = da.findGenotypeLabelTime(genotype, labels, tt, deltaTime) ;
                meanPIV_tt = snap.getMeanPIV() ;
                pivStack.vx(tidx,:,:) = meanPIV_tt.vx ;
                pivStack.vy(tidx,:,:) = meanPIV_tt.vy ;
                tidx = tidx + 1 ;
            end
            pivStack.npivx = size(meanPIV_tt.vy, 1) ;
            pivStack.npivy = size(meanPIV_tt.vy, 2) ;
            pivStack.x = Y0 ;
            pivStack.y = X0 ;
        end
        
        function t0V = getDynamicDataT0V(obj, options)
            % For dynamic embryo entry of the queriedSample
            method = 'pivlab' ;
            if isfield(options, 'pivmethod')
                method = options.pivmethod ;
            end
            
            error('add this function here')
        end
        
        function [XX, YY] = getPullbackPathlines(obj, pivStack, x0, y0, t0, options)
            % [XX, YY] = getPullbackPathlines(obj, pivStack, x0, y0, options)
            %
            % Parameters
            % ----------
            % obj : queriedSample class instance
            % pivStack : struct with two #TP x nX x nY float entries
            %   vx: #TP x nX x nY float
            %       velocities in x direction
            %   vy: #TP x nX x nY float
            %       velocities in y direction
            %   x : nX x nY numeric
            %       piv evaluation coordinates in x
            %   y : nX x nY numeric
            %       piv evaluation coordinates in y
            % x0 : n*m float array 
            %   x coordinates in pullback pixels to start pathlines at t0
            % y0 : n*m float array 
            %   y coordinates in pullback pixels to start pathlines at t0
            % t0 : 
            %   time at which to begin the pathlines, if member of
            %   options.timePoints. Otherwise, treated as an index into the PIV arrays
            % options : optional struct with fields
            %   struct passed to pullbackPathlines()
            % 
            % Example usage
            % -------------
            % 
            % 
            if nargin < 2
                error('Run qs.buildPIVStack() and pass pivStack to get pullback pathlines')
            end
            if nargin < 5
                options = struct() ;
            end
            if length(obj.meta.folders) == 1
                disp('Loading or computing pbPathlines for single embryo in queriedSample...')
                fnout = fullfile(obj.meta.folders{1}, 'pullbackPathlines.mat') ;
                if ~exist(fnout, 'file')
                    disp(['Computing pbPathlines & saving to disk: ' fnout])
                    [XX, YY] = pullbackPathlines(pivStack, x0, y0, t0, options) ;
                    XX = XX ./ pivStack.rescaleFactor ;
                    YY = YY ./ pivStack.rescaleFactor ;
                    save(fnout, 'XX', 'YY') ;
                else
                    disp(['Loading pbPathlines from disk: ' fnout])
                    load(fnout, 'XX', 'YY') ;
                end
            else
                disp('queriedSample contains multiple samples, computing pbPathlines without saving/loading...')
                [XX, YY] = pullbackPathlines(pivStack, x0, y0, t0, options) ;
                XX = XX ./ pivStack.rescaleFactor ;
                YY = YY ./ pivStack.rescaleFactor ;
            end
        end
        
        function pushForwardEmbryo(obj,XX,YY)
            
            if isempty(obj.pushForward)
                %         % our little mesh with PIV coordinates avs vertices
                %         mesh0 = struct() ;
                %         nU = size(X0v, 2) ;
                %         nV = size(Y0v, 1) ;
                %         X0t = X0v' ;
                %         Y0t = Y0v' ;
                %         uv = cat(3, X0t, Y0t) ;
                %         mesh0.f = defineFacesRectilinearGrid(uv, nU, nV) ;
                %         uu = reshape(uv(:, :, 1), [nU, nV]) ;
                %         vv = reshape(uv(:, :, 2), [nU, nV]) ;
                %         mesh0.u = [uu(:), vv(:)] ;
                %         mesh0.nU = nU ;
                %         mesh0.nV = nV ;
                %         % check it
                %         clf
                %         trisurf(triangulation(mesh0.f, mesh0.u(:, 1), mesh0.u(:, 2), 0*mesh0.u(:, 1)), 'facecolor', 'none')
                %         axis equal
                %     mesh0.v = zeros(size(mesh0.u, 1), 3) ;
                %     mesh0.v(:, 1) = xi(mesh0.u(:, 1), mesh0.u(:, 2)) ;
                %     mesh0.v(:, 2) = yi(mesh0.u(:, 1), mesh0.u(:, 2)) ;
                %     mesh0.v(:, 3) = zi(mesh0.u(:, 1), mesh0.u(:, 2)) ;

                % Larger mesh to use for interpolation
                mesh3d = read_ply_mod(fullfile(rootdir, 'embryo_geometry/embryo_rect_noglue.ply')) ;
                mesh3d.v = mesh3d.v * 0.2619 ; % convert to um

                mesh2d = read_ply_mod(fullfile(rootdir, 'embryo_geometry/rect_PIVImageScale.ply')) ;
                mesh3d.u = mesh2d.v ;

                assert(all(all(mesh2d.f == mesh3d.f)))

                % Project mesh0 into 3d using mesh3d
                xi = scatteredInterpolant(mesh3d.u(:, 2), mesh3d.u(:, 1), mesh3d.v(:, 1)) ;
                yi = scatteredInterpolant(mesh3d.u(:, 2), mesh3d.u(:, 1), mesh3d.v(:, 2)) ;
                zi = scatteredInterpolant(mesh3d.u(:, 2), mesh3d.u(:, 1), mesh3d.v(:, 3)) ;
                obj.pushForward.xi = xi ;
                obj.pushForward.yi = yi ;
                obj.pushForward.zi = zi ;
            end
        end
            
        
    end
    
end
