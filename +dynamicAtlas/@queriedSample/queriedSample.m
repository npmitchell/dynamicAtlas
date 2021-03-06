classdef queriedSample < handle
    %QUERIEDSAMPLE Object for manipulating subset of atlas
    %   Class for manipulating a queried subset of the dynamicAtlas
    %   Holds a subset of embryos in the atlas and offers methods for
    %   interrogating the properties of those embryos, such as loading the
    %   data, smoothed data, gradients, flow fields.
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
        piv = struct('vx', [], ...
            'vy', [], ...
            'rescaleFactor', 0.4)   % cell array of velocity fields
        meanPIV = struct('vx', [], ...
            'vy', [])           % cell array of average velocity field across sample
    end
    
    
    methods
        function obj = queriedSample(sampleStruct)
            obj.meta = sampleStruct ;
        end
        
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
                    if isfield(obj.meta, 'tiffpages')
                        disp(['Loading TIFF page ', ...
                            num2str(obj.meta.tiffpages(qq)), ...
                            ' from: ', tiff_fn])
                        try
                            dataCell{qq} = imread(tiff_fn, obj.meta.tiffpages(qq)) ;
                        catch
                            disp(['Could not load file. Incorrectly named? '  tiff_fn])
                            dataCell{qq} = imread(tiff_fn, obj.meta.tiffpages(qq)) ;
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
        
        function getPIV(obj)
            % Load PIV results for all entries in queriedSample collection
            nEmbryos = length(obj.meta.folders) ;
            rescaleFactor = obj.piv.rescaleFactor ;
            for qq = 1:nEmbryos
                if isa(obj.meta.times, 'cell')
                    % Load each in cell of timestamps
                    ntimestamps = obj.meta.nTimePoints(qq) ;
                    timestamps = 1:ntimestamps ;
                    vxcollection = zeros(46, 54, length(timestamps)) ;
                    vycollection = zeros(46, 54, length(timestamps)) ;
                    for pp = 1:(ntimestamps-1)
                        tiffpage = timestamps(pp) ;
                        try
                            pivfn = fullfile(obj.meta.folders{qq}, 'PIV_filtered', ...
                                sprintf('VeloT_medfilt_%06d.mat', tiffpage)) ;
                            tmp = load(pivfn) ;
                            vxcollection(:, :, pp) = tmp.VX / rescaleFactor ;
                            vycollection(:, :, pp) = tmp.VY / rescaleFactor ;
                        catch
                            disp('Warning: some timepoints have no velocity, returned these as NaN')
                        
                            vxcollection(:, :, pp) = NaN ;
                            vycollection(:, :, pp) = NaN ;
                        end
                    end
                    obj.piv.vx{qq} = vxcollection ;
                    obj.piv.vy{qq} = vycollection ;
                else
                    timestamp = obj.meta.tiffpages(qq) ;
                    pivfn = fullfile(obj.meta.folders{qq}, 'PIV_filtered', ...
                        sprintf('VeloT_medfilt_%06d.mat', timestamp)) ;
                    tmp = load(pivfn) ;
                    obj.piv.vx{qq} = tmp.VX / rescaleFactor ;
                    obj.piv.vy{qq} = tmp.VY / rescaleFactor ;
                end
            end
        end
        
        function meanPIV = getMeanPIV(obj)
            % Load PIV results for all entries in queriedSample collection
            nEmbryos = length(obj.meta.folders) ;
            vx = zeros(46, 54, nEmbryos) ;
            vy = zeros(46, 54, nEmbryos) ;
            rescaleFactor = obj.piv.rescaleFactor ;
            for qq = 1:nEmbryos
                if isa(obj.meta.times, 'cell')
                    error('handle this case of collections here -- see above in getPIV()')
                else
                    timestamp = obj.meta.tiffpages(qq) ;
                    pivfn = fullfile(obj.meta.folders{qq}, 'PIV_filtered', ...
                        sprintf('VeloT_medfilt_%06d.mat', timestamp)) ;
                    if exist(pivfn, 'file')
                        tmp = load(pivfn) ;
                        vx(:, :, qq) = tmp.VX / rescaleFactor ;
                        vy(:, :, qq) = tmp.VY / rescaleFactor ;
                    else
                        if timestamp == obj.meta.nTimePoints(qq)
                            disp('Skipping final timepoint from meanPIV')
                            vx(:, :, qq) = NaN ;
                            vy(:, :, qq) = NaN ;
                        else
                            disp(['PIV does not exist but tp is not final tp: ' pivfn])
                            vx(:, :, qq) = NaN ;
                            vy(:, :, qq) = NaN ;
                        end
                    end
                end
            end
            obj.meanPIV.vx = nanmean(vx, 3) ;
            obj.meanPIV.vy = nanmean(vy, 3) ;
            if nargout > 0
                meanPIV = obj.meanPIV ;
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
