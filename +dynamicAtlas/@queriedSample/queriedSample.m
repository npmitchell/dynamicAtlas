classdef queriedSample < handle
    %QUERIEDSAMPLE Object for manipulating subset of atlas
    %   Class for manipulating a queried subset of the dynamicAtlas
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
        meta        % struct with fields labels, folders, names, ...
                    %       embryoIDs, times, uncs, tiffpages, nTimePoints
        data        % cell array of images
        gradients = struct()
        smooth      % cell array of images
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
                        dataCell{qq} = imread(tiff_fn, obj.meta.tiffpages(qq)) ;
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
        
    end
    
end

