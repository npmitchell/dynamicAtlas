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
    % a = dynamicAtlas('./', {'WT'}) 
    % qs = a.findLabel('Eve')
    % qs.getData()
    % qs.getGradients('dx', 10)
    % qs.getSmooth()
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
                    if iscell(obj.meta.tiffpages)
                        thesePages = obj.meta.tiffpages{qq} ;
                    else
                        thesePages = obj.meta.tiffpages(qq) ;
                    end
                    if isfield(obj.meta, 'tiffpages') && ~iscell(obj.meta.tiffpages)
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
                    runMedianFilterOnPIVField(input_path, output_path, median_order);
                end
            end
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
            %   vx
            %   vy 
            %
            method = 'default';
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
                        vxcollection = zeros(86, 101, length(timestamps)) ;
                        vycollection = zeros(86, 101, length(timestamps)) ;
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
                                vxcollection(:, :, pp) = tmp.VX / (rescaleFactor *dt) ;
                                vycollection(:, :, pp) = tmp.VY / (rescaleFactor *dt) ;
                            else
                                pivfnRaw = fullfile(obj.meta.folders{qq}, 'PIVlab', ...
                                    sprintf('VeloT_fine_%06d.mat', tiffpage)) ;
                                pivfn = fullfile(obj.meta.folders{qq}, 'PIVlab_filtered', ...
                                    sprintf('VeloT_medfilt_fine_%06d.mat', tiffpage)) ;
                                tmp = load(pivfn) ;
                                if first
                                    load(pivfnRaw, 'convert_to_um_per_min', 'X0', 'Y0') ;
                                    obj.piv.X0 = X0 ;
                                    obj.piv.Y0 = Y0 ;
                                    first = false ;
                                else
                                    load(pivfnRaw, 'convert_to_um_per_min') ;
                                end
                                vxcollection(:, :, pp) = tmp.VX * convert_to_um_per_min ;
                                vycollection(:, :, pp) = tmp.VY * convert_to_um_per_min ;
                            end
                        catch
                            disp('Warning: some timepoints have no velocity, returned these as NaN')

                            vxcollection(:, :, pp) = NaN ;
                            vycollection(:, :, pp) = NaN ;
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
            % meanPIV : 
            %   mean flow across all samples in the queriedSample
            
            % default options
            method = 'default' ;
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
                nPages = nPages + length(obj.meta.tiffpages{qq}) ;
            end
            
            % preallocate the velocity information
            if strcmpi(method, 'default')
                vx = zeros(46, 54, nPages) ;
                vy = zeros(46, 54, nPages) ;
            else
                vx = zeros(86, 101, nPages) ;
                vy = zeros(86, 101, nPages) ;
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
                            vx(:, :, dmyk) = tmp.VX / (rescaleFactor * dt) ;
                            vy(:, :, dmyk) = tmp.VY / (rescaleFactor * dt) ;
                        else
                            if timestamp == obj.meta.nTimePoints(qq)
                                disp('Skipping final timepoint from meanPIV')
                                vx(:, :, dmyk) = NaN ;
                                vy(:, :, dmyk) = NaN ;
                            else
                                disp(['PIV does not exist but tp is not final tp: ' pivfn])
                                vx(:, :, dmyk) = NaN ;
                                vy(:, :, dmyk) = NaN ;
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
                            if strcmpi(method, 'pivlab_filtered') 
                                tmp = load(pivfn) ;
                                load(pivfnRaw, 'convert_to_um_per_min') ;
                                tmp.convert_to_um_per_min = convert_to_um_per_min ;
                            else
                                tmp = load(pivfn) ;
                            end
                            vx(:, :, dmyk) = tmp.VX * tmp.convert_to_um_per_min ;
                            vy(:, :, dmyk) = tmp.VY * tmp.convert_to_um_per_min ;
                            assert(any(~isnan(vx(:))))
                        else
                            if timestamp == obj.meta.nTimePoints(qq)
                                disp('Skipping final timepoint from meanPIV')
                                vx(:, :, dmyk) = NaN ;
                                vy(:, :, dmyk) = NaN ;
                            else
                                disp(['PIV does not exist but tp is not final tp: ' pivfn])
                                vx(:, :, dmyk) = NaN ;
                                vy(:, :, dmyk) = NaN ;
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
                obj.meanPIV.convert_to_um_per_min = tmp.convert_to_um_per_min ;
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
                obj.meanPIV.convert_to_um_per_min = (rescaleFactor * dt) ;
            end
            
            % Take mean
            obj.meanPIV.vx = mean(vx, 3, 'omitnan') ;
            obj.meanPIV.vy = mean(vy, 3, 'omitnan') ;
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
                    quiver(X0, Y0, squeeze(vx(:, :, dmyk)) * 10, ...
                        squeeze(vy(:, :, dmyk)) * 10, 0)
                    title(['page = ' num2str(dmyk)])
                    axis equal 
                    axis tight
                    pause(0.01)
                end
            end
            
        end
        
        function pivStack = buildPIVStack(obj, da, genotype, labels, deltaTime)
            % Build a stack of PIV flow fields over time
            
            if nargin < 5
                % half-width of search window for each timepoint
                deltaTime = 0.5 ;
            end
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
        
    end
    
end
