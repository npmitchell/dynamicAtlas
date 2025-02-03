function [ patchIm, imref, zeroID, MIP, SIP ] = ...
    texture_patch_to_image( FF, VV, TF, TV, IV, Options )
%TEXTURE_PATCH_TO_IMAGE This function generates an image from the texture
%mapping of a mesh triangulation with respect to an texture volume I. One 
%input is the "physical" or "real-space" mesh triangulation. This
%is the mesh triangulation you actually wish to display.  The other input
%is the "virtual" or "texture-space" mesh triangulation.  This input maps
%the physical mesh into texture space.  The virtual mesh must have the same
%number of facets as the physical mesh (each face needs a texture
%mapping!), but may have different connectivities (i.e. the physical mesh
%may be simply connected while the virtual mesh is multiply connected) or a
%different number of vertices (i.e. re-using the texturing from a single
%template face for multiple physical faces).
%
%   texture_patch_to_image( FF, VV, TF, TV, IV, Options )
%
%   Input Parameters:
%       - FF:       #Fx3 face connectivity list of the physical mesh
%       - VV:       #Vx2 vertex coordinate list of the physical mesh.
%                   NOTE: VV is provided in (X,Y) format
%       - TF:       #Fx3 face connectivity list of the virtual mesh
%       - TV:       #KxD texture coordinate list. Coordinates can be given
%                   as real pixel positions in the texture image or as a
%                   range [0..1]. D = (2,3). NOTE: TV is provided in
%                   (row, column)/(row, column, page) format
%       - IV:        The texture image volume. This function supports:
%                   Grayscale:  [I1 x I2] (2D)
%                   Grayscale:  [I1 x I2 x I3] (3D)
%                   RGB:    [I1 x I2 x 3] (2D)
%                   RGB:    { [I1 x I2 x I3] } x 3 cell array (3D)
%
%       - Options:  Structure containing the standard options for a
%                   textured surface patch, such as EdgeColor, EdgeAlpha,
%                   etc.  See MATLAB documentation for more information.
%
%       - Options.imSize:       The size of the output image
%       - Options.baseSize:     The side length in pixels of the smallest
%                               side of the output image when using the
%                               tight mesh bounding box
%       - Options.xLim:         The x-bounds of the output image
%       - Options.yLim:         The y-bounds of the output image
%       - Options.pixelSearch:  The method of searching for the faces
%                               containing pixel centers
%                                   - 'AABB' (requires GPToolBox)
%                                   - 'Default' (MATLAB built-ins, faster than AABB)
%       - Options.numLayers:    The number of onion layers to create
%                               Format is [ (num +), (num -) ]
%       - Options.layerSpacing: The spacing between adjacent onion layers
%                               in units of pixels
%       - Options.smoothIter:   Number of iterations of Laplacian mesh
%                               smoothing to run on the mesh prior to
%                               vertex normal displacement (requires
%                               GPToolBox) (Default is 0)
%       - Options.vertexNormal: User supplied vertex unit normals to the
%                               texture triangulation
%       - Options.Interpolant:  A pre-made texture image volume interpolant
%                               For RGB images the interpolant must be a
%                               three-element cell array with a separate
%                               interpolant for each color channel
%       - Options.isRGB:        A boolean indicating if the input image
%                               texture volume is RGB
%       - Options.isFalseColor: bool, whether to colorize with unique
%                               color for each channel
%       - Options.scaleData:    A boolean indicating if the ouput image
%                               data should be scaled or left with their
%                               raw values (true)
%       - Options.falseColors : #channels x 3 list of colors for each
%                               channel
%       - Options.Imax:         float, maximum value for the data interpolant 
%                               object above which we clip the intensity
%       - Options.Imin:         float, minimum value for the data interpolant 
%                               object below which we clip the intensity
%       - Options.extrapolationMethod: 
%         'nearest' | 'linear' | 'nearest' | 'next' | 'previous' | 
%         'pchip' | 'cubic' | 'spline' | 'makima' | 'none', 
%                               what extrapolation to use outside of the 
%                               domain of data interpolation values
%
%
%   Output Parameters:
%       - patchIm:  The output image stack
%       - imref:    The image spatial referencing object
%       - zeroID:   patchIm(:,:,zeroID) is data layer zero
%       - MIP:      The maximum intensity projection of the output stack
%       - SIP:      The summed intensity projection of the output stack
%
%   by Dillon Cislo 10/25/2019, NPMitchell 2020
%   NPMitchell added extrapolationMethod options
%   NPMitchell added functionality for false color multiple channels
%   Todo: add Imax, Imin for clipping the interpolated intensities
%    

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

% Verify input triangulations ---------------------------------------------

% Check that the number of faces is consistent
if ( size(FF,1) ~= size(TF,1) )
    error('texture_patch_to_image:inputs', ...
        [ 'Number of real-space faces must match', ...
        'the number of texture-space faces' ] );
end

% Check that the faces are triangles
if ( size(FF,2) ~= 3 ) || ( size(TF,2) ~= 3 )
    error('texture_patch_to_image:inputs', ...
        'All mesh faces must be elements of a triangulation' );
end

% Check the dimensions of the physical vertices
if( size(VV,2) ~= 2 )
    error('texture_patch_to_image:inputs', ...
        'Real-space vertex coordinates must be 2D');
end

% Check the dimensions of the texture vertices
if ~( (size(TV,2) == 3) || (size(TV,2) == 2) )
    error('texture_patch_to_image:inputs', ...
        'Texture-space coordinates are improperly sized');
end

% Re-scale texture vertices to true pixel positions if necessary
if max(TV(:)) < 2 
    for i = size(TV,2)
        TV(:,i) = (size(IV,i)-1) .* TV(:,i) + 1;
    end
end

% Create triangulation objects
realTri = triangulation( FF, VV );
textureTri = triangulation( TF, TV );

% Extract input options ---------------------------------------------------

% Use default options if none are supplied
if nargin < 6, Options = struct(); end

% Extract the output image x-bounds
if isfield( Options, 'xLim' )
    xLim = Options.xLim;
else
    xLim = [ min(VV(:,1)) max(VV(:,1)) ];
end

% Extract the output image y-bounds
if isfield( Options, 'yLim' )
    yLim = Options.yLim;
else
    yLim = [ min(VV(:,2)) max(VV(:,2)) ];
end

% Extract base image size
if isfield( Options, 'baseSize')
    baseSize = Options.baseSize;
else
    baseSize = 512;
end

% Extract the output image size
if isfield( Options, 'imSize' )
    imSize = Options.imSize;
else
    bndRatio = [ yLim(2)-yLim(1), xLim(2)-xLim(1) ];
    bndRatio = bndRatio ./ min(bndRatio);
    imSize = ceil( baseSize .* bndRatio );
end

% Construct the image spatial referencing object
if nargout > 1
    imref = imref2d( imSize, xLim, yLim );
end

% Determine the pixel center search method
if isfield( Options, 'pixelSearch' )
    pixelSearch = Options.pixelSearch;
    % Check that GPToolBox is on the path
    if strcmp(pixelSearch, 'AABB')
        if ~exist('in_element_aabb', 'file')
            warning('GPToolbox was not found. Using default settings.');
            pixelSearch = 'Default';
        end
    end
else
    pixelSearch = 'Default';
end

% Determine if any onion layers are to be produced
makePosLayers = false;
makeNegLayers = false;
makeOnion = false;
if isfield( Options, 'numLayers' )
    numLayers = Options.numLayers;
    if (numLayers(1) > 0)
        makePosLayers = true ;
    elseif numLayers(1) < 0
        error('number of layers in positive direction is negative')
    end
    if (abs(numLayers(2)) > 0), makeNegLayers = true; end
    if ( makePosLayers || makeNegLayers ), makeOnion = true; end
else
    numLayers = [0 0];
end

% Determine the location of the data layer zero
if nargout > 2
    zeroID = numLayers(2) + 1;
end

% Determine the onion layer spacing
if isfield( Options, 'layerSpacing' )
    layerSpacing = Options.layerSpacing;
    try
        assert(length(layerSpacing)==1)
    catch
        error('Options.layerSpacing must be an int/float value')
    end
else
    layerSpacing = 5;
end

% Determine how many iterations of Laplacian mesh smoothing to run on the
% input mesh prior to vertex displacement along normal vectors
if isfield( Options, 'smoothIter' )
    smoothIter = Options.smoothIter;
    % Check that GPToolBox is on the path
    if smoothIter > 0
        if ~exist('laplacian_smooth', 'file')
            warning('GPToolBox was not found. Using default settings.');
            smoothIter = 0;
            smoothMesh = false;
        else
            smoothMesh = true;
        end
    else
        smoothMesh = false;
    end
else
    smoothIter = 0;
    smoothMesh = false;
end

% Process texture vertex normals
if isfield( Options, 'vertexNormal' )
    vN = Options.vertexNormal;
else
    vN = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine if input is RGB or grayscale or falseColor
if isfield( Options, 'isRGB' )
    isRGB = Options.isRGB;
    Options = rmfield( Options, 'isRGB' );
    isFalseColor = false ;
elseif isfield(Options, 'isFalseColor')
    isFalseColor = Options.isFalseColor ;
    Options = rmfield( Options, 'isFalseColor' ) ;
    isRGB = false ;
else
    disp('determining if isFalseColor or isRGB')
    if ( iscell(IV) && length(IV) == 3 )
        % Assume that if IV is a length3 cell, its channels are RGB
        isRGB = true;
        isFalseColor = false ;
    elseif iscell(IV)
        isRGB = false ;
        isFalseColor = true ;
    else
        isRGB = false;
        isFalseColor = false ;
    end
end
if ~isRGB && isfield(Options, 'isFalseColor')
    isFalseColor = Options.isFalseColor ;
    Options = rmfield( Options, 'isFalseColor' ) ;
end

% Determine colors for falseColor option (when IV is a cell)
if isFalseColor
    if isfield(Options, 'falseColors')
        falseColors = Options.falseColors ;
        % Check that the correct number of colors is given
        if iscell(falseColors)
            if length(falseColors) ~= length(IV)
                error('Must pass same number of colors as number of color channels in falseColor rendering mode')
            end
        else
            if size(falseColors, 1) ~= length(IV)
                error('Must pass same number of colors as number of color channels in falseColor rendering mode')
            else
                % Convert falseColors to cell array if not done already
                if ~iscell(falseColors)
                    fC = cell(length(IV), 1) ;
                    for kk = 1:length(IV)
                        fC{kk} = falseColors(kk, :) ;
                    end
                end
                falseColors = fC ;
            end
        end
        Options = rmfield(Options, 'falseColors') ;
    else
        if length(IV) == 1
            disp('There is only one channel --> unpack it')
            IV = IV{1} ;
            isFalseColor = false ;
        elseif length(IV) == 2
            disp('falseColors not supplied --> Using default 2-color channels of red and cyan')
            falseColors{1} = [1, 0, 0] ;
            falseColors{2} = [0, 1, 1] ;
        elseif length(IV) == 3
            falseColors{1} = [1, 0, 0] ;
            falseColors{2} = [0, 1, 0] ;
            falseColors{3} = [0, 0, 1] ; 
        else
            error('For texture spaces with more than 3 channels, please supply colors for each channel')
        end
    end
end

% Determine intensity limits from Options
if isfield( Options, 'Imax' )
    Imax = Options.Imax;
    Options = rmfield( Options, 'Imax' );
else
    Imax = Inf ;
end
if isfield( Options, 'Imin' )
    Imin = Options.Imin;
    Options = rmfield( Options, 'Imin' );
else
    Imin = -Inf ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unpack extrapolationMethod
extrapolationMethod = 'nearest' ;
if isfield( Options, 'extrapolationMethod')
    extrapolationMethod = Options.extrapolationMethod ;
    Options = rmfield( Options, 'extrapolationMethod') ;
end











% Determine if output data should be scaled
if isfield( Options, 'scaleData' )
    scaleData = Options.scaleData;
else
    scaleData = true;
end

% Validate the input texture image volume ---------------------------------
if ~isempty(IV) % Allow users to supply pre-made interpolant by skipping here
    
    if isRGB
        
        if size(TV,2) == 2 % 2D
            
            if ~( (ndims(IV) == 3) && (size(IV,3) == 3) )
                error('texture_patch_to_image:inputs', ...
                    'Invalid texture image input');
            end
            
        else % 3D
            
            if ~( iscell(IV) && (numel(IV) == 3) )
                error('texture_patch_to_image:inputs', ...
                    'Invalid texture image input');
            else
                
               goodI = (ndims(IV{1}) == 3) && ...
                    isequal(size(IV{1}), size(IV{2})) && ...
                    isequal(size(IV{1}), size(IV{3})) && ...
                    isequal(size(IV{2}), size(IV{3}));
                
                if ~goodI
                    error('texture_patch_to_image:inputs', ...
                        'Invalid texture image input');
                end
                
            end
            
        end
        
    elseif isFalseColor
        
        if ~iscell(IV)
            error('texture_patch_3d:inputs', ...
                'Invalid texture image input: IV should be cell if Options.isFalseColor==true');
        end
        
    else
        
        if ~( (size(TV,2) == 3) && (ndims(IV) == 3) )
            if ~( (size(TV,2) == 2) && ismatrix(IV) )
                error('texture_patch_to_image:inputs', ...
                    'Invalid texture vertex/image input');
            end
        end
        
    end
    
end

% Create the interpolant object from the input image object ---------------
if isfield( Options, 'Interpolant' )
    
    IVI = Options.Interpolant;
    
    if isRGB
        
        % Extract each interpolant object for R, G, B separately
        if ~( iscell(IVI) && (numel(IVI) == 3) )
            error('texture_patch_to_image:inputs', ...
                'Invalid texture image interpolation object');
        end
        
        IIg = IVI{2}; IIb = IVI{3}; IVr = IVI{1};
        
    elseif isFalseColor
        
        % Extract each interpolant object for each face color
        if ~( iscell(IVI) && (numel(IVI) == length(IV)) )
            error('texture_patch_3d:inputs', ...
                'Invalid texture image interpolation object: must have same #elements as IV');
        end
        
        IVIfc = cell(size(IV)) ;
        for kk = 1:length(IV)
            IVIfc{kk} = IVI{kk} ;
        end
        clearvars IVI
        
    end
        
else
    
    if isRGB
        
        IVr = griddedInterpolant(single(IV{1}), 'cubic', extrapolationMethod);
        IIg = griddedInterpolant(single(IV{2}), 'cubic', extrapolationMethod);
        IIb = griddedInterpolant(single(IV{3}), 'cubic', extrapolationMethod);
                
    elseif isFalseColor
        % Each input texture channel receives its own interpolant
        IVIfc = cell(size(IV)) ;
        for kk = 1:length(IV)
            IVIfc{kk} = griddedInterpolant(single(IV{kk}), 'cubic', extrapolationMethod);
        end
        
    else
    
        IVI = griddedInterpolant(single(IV), 'cubic', extrapolationMethod);
        
    end
    
end

%--------------------------------------------------------------------------
% CREATE TEXTURE PATCH IMAGE FOR DATA LAYER ZERO
%--------------------------------------------------------------------------

% Calculate the centers of the pixels in real-space ----------------------

% The real-space spacing between adjacent columns
dx = ( xLim(2) - xLim(1) ) / imSize(2);

% The real-space spacing between adjacent rows
dy = ( yLim(2) - yLim(1) ) / imSize(1);

% The coordinates of the centers of pixels in real-space
[X, Y] = meshgrid( (xLim(1)+dx/2):dx:(xLim(2)-dx/2), ...
    (yLim(1)+dy/2):dy:(yLim(2)-dy/2) );

XY = [ X(:), Y(:) ];

% Find the faces containing the specified pixel centers in real-space and
% the barycentric interpolation coordinates of the pixel centers with
% respect to those containing faces ---------------------------------------

if strcmp( pixelSearch, 'AABB' )
    
    % NOTE: Points not located within any face return 0
    [ pixFaces, ~, ~, ~ ] = in_element_aabb( VV, FF, XY, [], [], [] );
    pixFaces( pixFaces == 0 ) = NaN;
    
    pixBaryReal = nan( size(XY,1), 3 );
    pixBaryReal(~isnan(pixFaces), :) = cartesianToBarycentric( realTri, ...
        pixFaces(~isnan(pixFaces)), XY(~isnan(pixFaces), :) );

else
    
    % NOTE: Points not located within any face return NaN
    [ pixFaces, pixBaryReal ] = pointLocation( realTri, XY );
    
end

% Find the texture coordinates of the pixel centers -----------------------
pixTexture = nan( size(XY,1), size(TV,2) );
pixTexture(~isnan(pixFaces), :) = barycentricToCartesian( textureTri, ...
    pixFaces(~isnan(pixFaces)), pixBaryReal(~isnan(pixFaces), :) );

% Interpolate to find the texture intensity at the pixel center -----------

if isRGB % RGB
    
    PIR = zeros( size(XY,1), 1 );
    PIG = PIR; PIB = PIR;
    
    if (size(TV,2) == 2) % 2D
        
        % Red channel
        PIR(~isnan(pixFaces)) = IVr( ...
            pixTexture(~isnan(pixFaces), 1), ...
            pixTexture(~isnan(pixFaces), 2) );
        
        % Green channel
        PIG(~isnan(pixFaces)) = IIg( ...
            pixTexture(~isnan(pixFaces), 1), ...
            pixTexture(~isnan(pixFaces), 2) );
        
        % Blue channel
        PIB(~isnan(pixFaces)) = IIb( ...
            pixTexture(~isnan(pixFaces), 1), ...
            pixTexture(~isnan(pixFaces), 2) );
        
    else % 3D
        
        % Red channel
        PIR(~isnan(pixFaces)) = IVr( ...
            pixTexture(~isnan(pixFaces), 1), ...
            pixTexture(~isnan(pixFaces), 2), ...
            pixTexture(~isnan(pixFaces), 3) );
        
        % Green channel
        PIG(~isnan(pixFaces)) = IIg( ...
            pixTexture(~isnan(pixFaces), 1), ...
            pixTexture(~isnan(pixFaces), 2), ...
            pixTexture(~isnan(pixFaces), 3) );
        
        % Blue channel
        PIB(~isnan(pixFaces)) = IIb( ...
            pixTexture(~isnan(pixFaces), 1), ...
            pixTexture(~isnan(pixFaces), 2), ...
            pixTexture(~isnan(pixFaces), 3) );
        
    end
    
    % Resize to desired output dimensions
    PIR = reshape( PIR, imSize );
    PIG = reshape( PIG, imSize );
    PIB = reshape( PIB, imSize );
    
    % Construct RGB image
    patchIm = cat( 3, PIR, PIG, PIB );

elseif isFalseColor
    
    PIR = zeros( size(XY,1), 1 );
    PIG = PIR; PIB = PIR;
    
    
    if (size(TV,2) == 2) % 2D
        
        for kk = 1:length(IVIfc)
            i2add = IVIfc{kk}( ...
                pixTexture(~isnan(pixFaces), 1), ...
                pixTexture(~isnan(pixFaces), 2)) ;
            % Red channel
            PIR(~isnan(pixFaces)) = PIR(~isnan(pixFaces)) + ...
                i2add * falseColors{kk}(1) ;
            % Green channel
            PIG(~isnan(pixFaces)) = PIG(~isnan(pixFaces)) + ...
                i2add * falseColors{kk}(2) ;
            % Blue channel
            PIB(~isnan(pixFaces)) = PIB(~isnan(pixFaces)) + ...
                i2add * falseColors{kk}(3) ;
        end
        
    else % 3D
        
        for kk = 1:length(IVIfc)
            i2add = IVIfc{kk}( ...
                pixTexture(~isnan(pixFaces), 1), ...
                pixTexture(~isnan(pixFaces), 2), ...
                pixTexture(~isnan(pixFaces), 3) ) ;
            % Red channel
            PIR(~isnan(pixFaces)) = PIR(~isnan(pixFaces)) + ...
                i2add * falseColors{kk}(1) ;
            % Green channel
            PIG(~isnan(pixFaces)) = PIG(~isnan(pixFaces)) + ...
                i2add * falseColors{kk}(2) ;
            % Blue channel
            PIB(~isnan(pixFaces)) = PIB(~isnan(pixFaces)) + ...
                i2add * falseColors{kk}(3) ;
        end
        
    end    

    % Resize to desired output dimensions
    PIR = reshape( PIR, imSize );
    PIG = reshape( PIG, imSize );
    PIB = reshape( PIB, imSize );
    
    % Construct RGB image
    patchIm = cat( 3, PIR, PIG, PIB );
    
else % Grayscale
    
    patchIm = zeros( size(XY,1), 1 );
    
    if (size(TV,2) == 2) % 2D
        
        patchIm(~isnan(pixFaces)) = IVI( ...
            pixTexture(~isnan(pixFaces), 1), ...
            pixTexture(~isnan(pixFaces), 2) );
        
    else % 3D
        
        patchIm(~isnan(pixFaces)) = IVI( ...
            pixTexture(~isnan(pixFaces), 1), ...
            pixTexture(~isnan(pixFaces), 2), ...
            pixTexture(~isnan(pixFaces), 3 ));
        
    end
    
    % Resize to desired output dimensions
    patchIm = reshape( patchIm, imSize );
    
end

%--------------------------------------------------------------------------
% GENERATE ONION LAYERS
%--------------------------------------------------------------------------

if makeOnion
    
    % Smooth the texture triangulation if desired
    if smoothMesh
        % Hold the positions of the boundary vertices fixed under smoothing
        % to avoid mesh collapse
        bdyIDx = unique(realTri.freeBoundary);
        smoothV = laplacian_smooth( TV, TF, 'cotan', ...
            bdyIDx, 0.1, 'implicit', TV, smoothIter );
    else
        smoothV = textureTri.Points;
    end
    
    % Calculate the vertex unit normals of the texture triangulation
    if isempty(vN)
        % Default is angle weighted vertex normals

        % Note the minus sign - it is necessary to preserve a notion of
        % 'outward pointing' for the (row, column, page) format of the
        % texture vertices

        try

            mesh = struct();
            mesh.f = TF;
            mesh.v = smoothV;
            mesh.vn = -per_vertex_normals(smoothV, TF, 'Weighting', 'angle');
            vN = average_normals_with_neighboring_vertices(mesh, 0.5);
        catch

            % Alternative is uniform weighting of faces incident to each
            % vertex.
            vN = -vertexNormal( triangulation( TF, smoothV ) );

        end
        
    end
    
    %----------------------------------------------------------------------
    % Generate negative layers
    %----------------------------------------------------------------------
    
    if makeNegLayers
        
        % The stack holding the negative layers
        if isRGB || isFalseColor
            mStack = zeros( [ imSize 3 abs(numLayers(2)) ] );
        else
            mStack = zeros( [ imSize abs(numLayers(2)) ] );
        end
        
        for i = 1:abs(numLayers(2))
            
            % Determine the texture coordinates for the current layer
            OV = smoothV - i .* layerSpacing .* vN;
            
            % Construct a triangulation representation
            smoothTri = triangulation( TF, OV );
            
            % Find the texture coordinates of the pixel centers
            pixTexture = nan( size(XY,1), size(TV,2) );
            pixTexture(~isnan(pixFaces), :) = barycentricToCartesian( ...
                smoothTri, pixFaces(~isnan(pixFaces)), ...
                pixBaryReal(~isnan(pixFaces), :) );
            
            % Interpolate to find the texture intensity at the pixel center
            
            if isRGB % Color
                
                LIR = zeros( size(XY,1), 1 );
                LIG = LIR; LIB = LIR;
                
                if (size(TV,2) == 2) % 2D
                    
                    % Red channel
                    LIR(~isnan(pixFaces)) = IVI( ...
                        pixTexture(~isnan(pixFaces), 1), ...
                        pixTexture(~isnan(pixFaces), 2) );
                    
                    % Green channel
                    LIG(~isnan(pixFaces)) = IIg( ...
                        pixTexture(~isnan(pixFaces), 1), ...
                        pixTexture(~isnan(pixFaces), 2) );
                    
                    % Blue channel
                    LIB(~isnan(pixFaces)) = IIb( ...
                        pixTexture(~isnan(pixFaces), 1), ...
                        pixTexture(~isnan(pixFaces), 2) );
                    
                else % 3D
                    
                    % Red channel
                    LIR(~isnan(pixFaces)) = IVI( ...
                        pixTexture(~isnan(pixFaces), 1), ...
                        pixTexture(~isnan(pixFaces), 2), ...
                        pixTexture(~isnan(pixFaces), 3) );
                    
                    % Green channel
                    LIG(~isnan(pixFaces)) = IIg( ...
                        pixTexture(~isnan(pixFaces), 1), ...
                        pixTexture(~isnan(pixFaces), 2), ...
                        pixTexture(~isnan(pixFaces), 3) );
                    
                    % Blue channel
                    LIB(~isnan(pixFaces)) = IIb( ...
                        pixTexture(~isnan(pixFaces), 1), ...
                        pixTexture(~isnan(pixFaces), 2), ...
                        pixTexture(~isnan(pixFaces), 3) );
                    
                end
                
                % Resize to desired output dimensions
                LIR = reshape( LIR, imSize );
                LIG = reshape( LIG, imSize );
                LIB = reshape( LIB, imSize );
                                                
                % Add the current layer to the image stack
                mStack(:,:,:,i) = cat( 3, LIR, LIG, LIB );

            elseif isFalseColor

                LIR = zeros( size(XY,1), 1 );
                LIG = LIR; LIB = LIR;


                if (size(TV,2) == 2) % 2D

                    for kk = 1:length(IVIfc)
                        i2add = IVIfc{kk}( ...
                            pixTexture(~isnan(pixFaces), 1), ...
                            pixTexture(~isnan(pixFaces), 2)) ;
                        % Red channel
                        LIR(~isnan(pixFaces)) = LIR(~isnan(pixFaces)) + ...
                            i2add * falseColors(kk, 1) ;
                        % Green channel
                        LIG(~isnan(pixFaces)) = LIG(~isnan(pixFaces)) + ...
                            i2add * falseColors(kk, 2) ;
                        % Blue channel
                        LIB(~isnan(pixFaces)) = LIB(~isnan(pixFaces)) + ...
                            i2add * falseColors(kk, 3) ;
                    end

                else % 3D

                    for kk = 1:length(IVIfc)
                        i2add = IVIfc{kk}( ...
                            pixTexture(~isnan(pixFaces), 1), ...
                            pixTexture(~isnan(pixFaces), 2), ...
                            pixTexture(~isnan(pixFaces), 3) ) ;
                        % Red channel
                        LIR(~isnan(pixFaces)) = LIR(~isnan(pixFaces)) + ...
                            i2add * falseColors{kk}(1) ;
                        % Green channel
                        LIG(~isnan(pixFaces)) = LIG(~isnan(pixFaces)) + ...
                            i2add * falseColors{kk}(2) ;
                        % Blue channel
                        LIB(~isnan(pixFaces)) = LIB(~isnan(pixFaces)) + ...
                            i2add * falseColors{kk}(3) ;
                    end

                end    

                % Resize to desired output dimensions
                LIR = reshape( LIR, imSize );
                LIG = reshape( LIG, imSize );
                LIB = reshape( LIB, imSize );
                
                % Construct RGB image
                patchIm = cat( 3, LIR, LIG, LIB );

            else % Grayscale
                
                layerIm = zeros( size(XY,1), 1 );
                
                if size(TV,2) == 2 % 2D
                    
                    layerIm(~isnan(pixFaces)) = IVI( ...
                        pixTexture(~isnan(pixFaces), 1), ...
                        pixTexture(~isnan(pixFaces), 2) );
                    
                else % 3D
                    
                    layerIm(~isnan(pixFaces)) = IVI( ...
                        pixTexture(~isnan(pixFaces), 1), ...
                        pixTexture(~isnan(pixFaces), 2), ...
                        pixTexture(~isnan(pixFaces), 3) );
                    
                end
                
                % Resize to desired output dimensions
                layerIm = reshape( layerIm, imSize );
                
                % Add the current layer image to the stack
                mStack(:,:,i) = layerIm;
                
            end
            
        end
        
        if isRGB || isFalseColor
            
            % Flip stack to reflect proper spatial ordering
            mStack = flip( mStack, 4 );
            
            % Concatenate the negative layers and data layer zero
            patchIm = cat( 4, mStack, patchIm );
            
        else
            
            % Flip stack to reflect proper spatial ordering
            mStack = flip( mStack, 3 );
            
            % Concatenate the negative layers and data layer zero
            patchIm = cat( 3, mStack, patchIm );
            
        end
        
    end
    
    if makePosLayers
        
        % The stack holding the positive layers
        if isRGB || isFalseColor
            pStack = zeros( [ imSize 3 abs(numLayers(1)) ] );
        else
            pStack = zeros( [ imSize abs(numLayers(1)) ] );
        end
        
        for i = 1:abs(numLayers(1))
            
            % Determine the texture coordinates for the current layer
            OV = smoothV + i .* layerSpacing .* vN;
            
            % Construct a triangulation representation
            smoothTri = triangulation( TF, OV );
            
            % Find the texture coordinates of the pixel centers
            pixTexture = nan( size(XY,1), size(TV,2) );
            pixTexture(~isnan(pixFaces), :) = barycentricToCartesian( ...
                smoothTri, pixFaces(~isnan(pixFaces)), ...
                pixBaryReal(~isnan(pixFaces), :) );
            
            % Interpolate to find the texture intensity at the pixel center
            
            if isRGB % Color
                
                LIR = zeros( size(XY,1), 1 );
                LIG = LIR; LIB = LIR;
                
                if (size(TV,2) == 2) % 2D
                    
                    % Red channel
                    LIR(~isnan(pixFaces)) = IVI( ...
                        pixTexture(~isnan(pixFaces), 1), ...
                        pixTexture(~isnan(pixFaces), 2) );
                    
                    % Green channel
                    LIG(~isnan(pixFaces)) = IIg( ...
                        pixTexture(~isnan(pixFaces), 1), ...
                        pixTexture(~isnan(pixFaces), 2) );
                    
                    % Blue channel
                    LIB(~isnan(pixFaces)) = IIb( ...
                        pixTexture(~isnan(pixFaces), 1), ...
                        pixTexture(~isnan(pixFaces), 2) );
                    
                else % 3D
                    
                    % Red channel
                    LIR(~isnan(pixFaces)) = IVI( ...
                        pixTexture(~isnan(pixFaces), 1), ...
                        pixTexture(~isnan(pixFaces), 2), ...
                        pixTexture(~isnan(pixFaces), 3) );
                    
                    % Green channel
                    LIG(~isnan(pixFaces)) = IIg( ...
                        pixTexture(~isnan(pixFaces), 1), ...
                        pixTexture(~isnan(pixFaces), 2), ...
                        pixTexture(~isnan(pixFaces), 3) );
                    
                    % Blue channel
                    LIB(~isnan(pixFaces)) = IIb( ...
                        pixTexture(~isnan(pixFaces), 1), ...
                        pixTexture(~isnan(pixFaces), 2), ...
                        pixTexture(~isnan(pixFaces), 3) );
                    
                end
                
                % Resize to desired output dimensions
                LIR = reshape( LIR, imSize );
                LIG = reshape( LIG, imSize );
                LIB = reshape( LIB, imSize );
                
                % Add the current layer to the image stack
                pStack(:,:,:,i) = cat( 3, LIR, LIG, LIB );
                
            elseif isFalseColor

                LIR = zeros( size(XY,1), 1 );
                LIG = LIR; LIB = LIR;


                if (size(TV,2) == 2) % 2D

                    for kk = 1:length(IVIfc)
                        i2add = IVIfc{kk}( ...
                            pixTexture(~isnan(pixFaces), 1), ...
                            pixTexture(~isnan(pixFaces), 2)) ;
                        % Red channel
                        LIR(~isnan(pixFaces)) = LIR(~isnan(pixFaces)) + ...
                            i2add * falseColors(kk, 1) ;
                        % Green channel
                        LIG(~isnan(pixFaces)) = LIG(~isnan(pixFaces)) + ...
                            i2add * falseColors(kk, 2) ;
                        % Blue channel
                        LIB(~isnan(pixFaces)) = LIB(~isnan(pixFaces)) + ...
                            i2add * falseColors(kk, 3) ;
                    end

                else % 3D

                    for kk = 1:length(IVIfc)
                        i2add = IVIfc{kk}( ...
                            pixTexture(~isnan(pixFaces), 1), ...
                            pixTexture(~isnan(pixFaces), 2), ...
                            pixTexture(~isnan(pixFaces), 3) ) ;
                        % Red channel
                        LIR(~isnan(pixFaces)) = LIR(~isnan(pixFaces)) + ...
                            i2add * falseColors{kk}(1) ;
                        % Green channel
                        LIG(~isnan(pixFaces)) = LIG(~isnan(pixFaces)) + ...
                            i2add * falseColors{kk}(2) ;
                        % Blue channel
                        LIB(~isnan(pixFaces)) = LIB(~isnan(pixFaces)) + ...
                            i2add * falseColors{kk}(3) ;
                    end

                end    

                % Resize to desired output dimensions
                LIR = reshape( LIR, imSize );
                LIG = reshape( LIG, imSize );
                LIB = reshape( LIB, imSize );
                
                % Construct RGB image
                patchIm = cat( 3, LIR, LIG, LIB );
                
                % Add the current layer to the image stack
                pStack(:,:,:,i) = cat( 3, LIR, LIG, LIB );

            else % Grayscale
                
                layerIm = zeros( size(XY,1), 1 );
                
                if size(TV,2) == 2
                    
                    layerIm(~isnan(pixFaces)) = IVI( ...
                        pixTexture(~isnan(pixFaces), 1), ...
                        pixTexture(~isnan(pixFaces), 2) );
                    
                else
                    
                    layerIm(~isnan(pixFaces)) = IVI( ...
                        pixTexture(~isnan(pixFaces), 1), ...
                        pixTexture(~isnan(pixFaces), 2), ...
                        pixTexture(~isnan(pixFaces), 3) );
                    
                end
                
                % Resize to desired output dimensions
                layerIm = reshape( layerIm, imSize );
                
                % Add the current layer image to the stack
                pStack(:,:,i) = layerIm;
                
            end
            
        end
        
        % Concatenate the positive layers and data layer zero
        if isRGB || isFalseColor
            patchIm = cat( 4, patchIm, pStack );
        else
        	patchIm = cat( 3, patchIm, pStack );
        end
        
    end
    
end

%--------------------------------------------------------------------------
% FORMAT OUTPUT
%--------------------------------------------------------------------------

% Calculate output stack MIP ----------------------------------------------
if nargin > 3
    if isRGB || isFalseColor
        if scaleData
            MIP = cat( 3, ...
                mat2gray( max( squeeze(patchIm(:,:,1,:)), [], 3 ) ), ...
                mat2gray( max( squeeze(patchIm(:,:,2,:)), [], 3 ) ), ...
                mat2gray( max( squeeze(patchIm(:,:,3,:)), [], 3 ) ) );
        else
            MIP = cat( 3, ...
                max( squeeze(patchIm(:,:,1,:)), [], 3 ), ...
                max( squeeze(patchIm(:,:,2,:)), [], 3 ), ...
                max( squeeze(patchIm(:,:,3,:)), [], 3 ) );
        end
    else
        if scaleData
            MIP = mat2gray(max( patchIm, [], 3 ));
        else
            MIP = max( patchIm, [], 3 );
        end
    end
end

% Calculate output stack SIP ----------------------------------------------
if nargin > 4
    if isRGB
        if scaleData
            SIP = cat( 3, ...
                mat2gray( sum( squeeze(patchIm(:,:,1,:)), 3 ) ), ...
                mat2gray( sum( squeeze(patchIm(:,:,2,:)), 3 ) ), ...
                mat2gray( sum( squeeze(patchIm(:,:,3,:)), 3 ) ) );
        else
            SIP = cat( 3, ...
                sum( squeeze(patchIm(:,:,1,:)), 3 ), ...
                sum( squeeze(patchIm(:,:,2,:)), 3 ), ...
                sum( squeeze(patchIm(:,:,3,:)), 3 ) );
        end
    else
        if scaleData
            SIP = mat2gray(sum( patchIm, 3 ));
        else
            SIP = sum( patchIm, 3 );
        end
    end
end

% Re-scale output image stack ---------------------------------------------

if scaleData
    
    % Re-scale each channel separately
    if isRGB
        
        for c = 1:3
            
            PIC = patchIm(:,:,c,:);
            
            limits = [ min(PIC(:)) max(PIC(:)) ];
            
            if limits(2) ~= limits(1)
                delta = 1 / (limits(2)-limits(1));
                patchIm(:,:,c,:) = delta .* ( PIC - limits(1) );
            end
            
        end
        
    else
        
        limits = [ min(patchIm(:)) max(patchIm(:)) ];
        
        if limits(2) ~= limits(1)
            delta = 1 / (limits(2) - limits(1));
            patchIm = delta .* ( patchIm - limits(1) );
        end
        
    end
    
    % Make sure all values are on the range [0, 1]
    patchIm( patchIm < 0 ) = 0;
    patchIm( patchIm > 1 ) = 1;
    
end


end

