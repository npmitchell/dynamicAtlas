function texture_patch_2d( FF, VV, TF, TV, I, Options)
%TEXTURE_PATCH_2D This function will display a triangulated mesh, like the
%MATLAB function "patch", with a texture mapping given by the 2D image I.
%Based in part on the function "patcht" by D. Kroon. One set of inputs to
%this function is the "physical" or "real-space" mesh triangulation. This
%is the mesh triangulation you actually wish to display.  The other input
%is the "virtual" or "texture-space" mesh triangulation.  This input maps
%the physical mesh into texture space.  The virtual mesh must have the same
%number of facets as the physical mesh (each face needs a texture
%mapping!), but may have different connectivities (i.e. the physical mesh
%may be simply connected while the virtual mesh is multiply connected) or a
%different number of vertices (i.e. re-using the texturing from a single
%template face for multiple physical faces).
%
%   texture_patch_2d( FF, VV, TF, TV, I, Options );
%
%   Input Parameters:
%       - FF:       #Fx3 face connectivity list of the physical mesh
%       - VV:       #Vx3 vertex coordinate list of the physical mesh. NOTE:
%                   VV is provided in (X,Y,Z) format!
%       - TF:       #Fx3 face connectivity list of the virtual mesh
%       - TV:       #Kx2  texture coordinate list. Coordinates can be given
%                   as real pixel positions in the texture image or as a
%                   range [0..1]. NOTE: TV is provided in (row, column)
%                   format!
%       - I:        The 2D texture image.  RGB: [I1 x I2 x 3] or Grayscale
%                   [I1 x I2]
%       - Options:  Structure containing the standard options for a
%                   textured surface patch, such as EdgeColor, EdgeAlpha,
%                   etc.  See MATLAB documentation for more information.
%
%       - Options.PSize:	Special option, defines the image texture size
%       for each individual polygon.  A lower number gives a more 
%       pixelated texture {64}
%
%       - Options.ApplyAmbientOcclusion:    Determines if the texture
%       colors should be modified by the ambient occlusion of the 
%       underlying triangulation {'false'}
%
%       - Options.AmbientOcclusion:    #Vx1 list of ambient occlusion
%       values {[]}
%
%       - Options.AmbientOcculsionFactor:   A scalar between [0,1]
%       determining how strongly the ambient occlusion modifies the 
%       texture mapping {1}
%
%       - Options.AmbientOcculsionSamples:  The number of samples to use
%       when computing the ambient occlusion {1000}
%
%       - Options.Unoriented:   Treat the surface as unoriented when
%       applying ambient occlusion {'false'}
%
%       - Options.AddLights:    Will add lighting to the plot {'false'}
%
%
%   by Dillon Cislo 08/13/2019


%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

% Patch FaceColor MUST be a texture.
Options.FaceColor = 'texturemap';

% Check that the number of faces is consistent
if ( size(FF,2) ~= size(TF,2) )
    error('texture_patch_2d:inputs', ...
        [ 'Number of real-space faces must match', ...
        'the number of texture-space faces' ] );
end

% Check the dimensions of the physical vertices
if ( size(VV,2) ~= 3 )
    if( size(VV,2) == 2 )
        VV = [ VV zeros( size(VV,1), 1 ) ];
    else
        error('texture_patch_2d:inputs', ...
            'Invalid real-space vertex input');
    end
end

% Check the dimensions of the texture vertices
if ( size(TV,2) ~= 2 )
    error('texture_patch_2d:inputs', 'Invalid texture vertex input');
end

% Check the dimensions of the texure image
if ( (~ismatrix(I)) && (ndims(I) ~= 3) )
    error('texture_patch_2d:inputs', 'Invalid texture image input');
end

% Detect if the texture image is grayscale or RGB
switch size(I,3)
    case 1
        iscolor = false;
    case 3
        iscolor = true;
    otherwise
        error('texture_patch_2d:inputs', 'Invalid texture image input');
end

% Re-scale texture vertices to true pixel positions if necessary
if ( max(TV(:)) < 2 )
    TV = [ (size(I,1)-1) .* TV(:,1) + 1, (size(I,2)-1) .* TV(:,2) + 1 ];
end

% Split texture image into separate R,G,B matrices to allow for fast linear
% indexing
Ir = I(:,:,1);
if iscolor, Ig = I(:,:,2); Ib = I(:,:,3); end

% Validate Optional Inputs ------------------------------------------------

% Size of the texture image used for each triangle
if isfield( Options, 'PSize' )
    sizep = round(Options.PSize(1));
    Options = rmfield(Options, 'PSize');
else
    sizep = 64;
end

if isfield( Options, 'ApplyAmbientOcclusion' )
    applyAO = Options.ApplyAmbientOcclusion;
    Options =rmfield(Options, 'ApplyAmbientOcclusion');
else
    applyAO = false;
end

if isfield( Options, 'AmbientOcclusion' )
    AO = Options.AmbientOcclusion;
    Options = rmfield( Options, 'AmbientOcclusion' );
    if ~isempty(AO), applyAO = true; end
else
    AO = [];
end

if isfield( Options, 'AmbientOcclusionFactor' )
    AOFactor = Options.AmbientOcclusionFactor;
    Options = rmfield( Options, 'AmbientOcclusionFactor' );
else
    AOFactor = 1;
end

if isfield( Options, 'AmbientOcclusionSamples' )
    AOSamples = Options.AmbientOcclusionSamples;
    Options = rmfield( Options, 'AmbientOcclusionSamples' );
else
    AOSamples = 1000;
end

if isfield( Options, 'Unoriented' )
    unoriented = Options.Unoriented;
    Options = rmfield( Options, 'Unoriented' );
else
    unoriented = false;
end

if isfield( Options, 'AddLights' )
    addLights = Options.AddLights;
    Options = rmfield( Options, 'AddLights' );
else
    addLights = false;
end

%--------------------------------------------------------------------------
% AMBIENT OCCLUSION HANDLING
%--------------------------------------------------------------------------

if applyAO
    
    if isempty(AO)
        
        % Calculate physical mesh vertex normals
        % NOTE: MATLAB uses backwards normals for graphics...
        VN = per_vertex_normals( VV, FF, 'Weighting', 'angle' );
        
        % Calculate ambient occulsion
        AO = ambient_occlusion(VV, FF, VV, -VN, AOSamples);
        if unoriented
            AO = min( AO, ambient_occlusion(VV, FF, VV, VN, AOSamples) );
        end
        
    end
    
end

%--------------------------------------------------------------------------
% FIND PATCH INTERPOLATION VALUES
%--------------------------------------------------------------------------

% The patch used for every triangle
% For simplicity these will be the same for each triangle regardless of
% size in either the physical or texture space. A smarter algorithm would
% probably take these features into account
Jr = zeros( (sizep+1), (sizep+1), class(I) );
if iscolor, Jg = Jr; Jb = Jr; end

% Linear indices of the 2D image associated with the triangle patch
jind = (sizep+1)^2:-1:1;

% Barycentric interpolation values for query points in the triangle patch
[ lambda1, lambda2, lambda3 ] = ...
    calculateBarycentricInterpolationValues( sizep );

%--------------------------------------------------------------------------
% CREATE THE SURFACE PLOT
%--------------------------------------------------------------------------

hold on

% Loop through all of the faces of the mesh
for i = 1:size(FF,1)
    
    % The currrent face's physical vertex coordinates
    V = VV( FF(i,:), : );
    
    % The current face's texture vertex coordinates
    tV = TV( TF(i,:), : );
    
    % Define the physical face in 'surface' format
    x=[V(1,1) V(2,1); V(3,1) V(3,1)];
    y=[V(1,2) V(2,2); V(3,2) V(3,2)];
    z=[V(1,3) V(2,3); V(3,3) V(3,3)];
    
    % Define the texture coordinates of the 'surface'
    xy = [ tV(1,:); tV(2,:); tV(3,:); tV(3,:) ];
    
    % Calculate the texture interpolation coordinates ---------------------
    
    pos(:,1) = xy(1,1)*lambda1 + xy(2,1)*lambda2 + xy(3,1)*lambda3;
    pos(:,2) = xy(1,2)*lambda1 + xy(2,2)*lambda2 + xy(3,2)*lambda3;
    
    % Round interpolation coordinates to true pixel positions
    pos = round(pos);
    pos = max(pos,1);
    pos(:,1) = min( pos(:,1), size(I,1) );
    pos(:,2) = min( pos(:,2), size(I,2) );
    
    % Convert coordinates to linear indices
    posind = (pos(:,1)-1) + (pos(:,2)-1) * size(I,1) + 1;
    
    % Map texture to surface image ----------------------------------------
    Jr(jind) = Ir(posind);
    J(:,:,1) = Jr; 
    if(iscolor)
        Jg(jind) = Ig(posind); 
        Jb(jind) = Ib(posind);
        J(:,:,2) = Jg; 
        J(:,:,3) = Jb;
    end
    
    % Apply ambient occlusion ---------------------------------------------
    if applyAO
        
        % Interpolate the ambient occlusion values on vertices
        AOMat = AO(FF(i,1))*lambda1 + AO(FF(i,2))*lambda2 + ...
            AO(FF(i,3))*lambda3;
        
        % Reshape to the size of the surface image
        AOMat = reshape(AOMat, sizep+1, sizep+1);
        
        % Change data type if necessary
        if isa( I, 'uint8' )
            J = double(J);
            Jr = double(Jr) ./ 255;
            if iscolor
                Jg = double(Jg) ./ 255;
                Jb = double(Jb) ./ 255;
            end
        end
        
        Jr = (1-AOFactor) .* Jr + AOFactor .* (Jr .* (1-AOMat));
        J(:,:,1) = Jr;
        if iscolor
            Jg = (1-AOFactor) .* Jg + AOFactor .* (Jg .* (1-AOMat));
            Jb = (1-AOFactor) .* Jb + AOFactor .* (Jb .* (1-AOMat));
            J(:,:,2) = Jg;
            J(:,:,3) = Jb;
        end
        
    end
    
    % Show the surface ----------------------------------------------------
    surface( x, y, z, J, Options );
    
end

hold off

%--------------------------------------------------------------------------
% Lighting Handling
%--------------------------------------------------------------------------

if addLights
    
    % Extract the center of the current axis bounding box
    XLim = get(gca, 'XLim');
    YLim = get(gca, 'YLim');
    ZLim = get(gca, 'ZLim');
    
    cen = [ mean(XLim) mean(YLim) mean(ZLim) ];
    
    light( 'Position', [cen(1) cen(2) 10*(max(ZLim)-min(ZLim))+cen(3)], ...
        'Style', 'local', 'Color', [1 1 1]/3 );
    light( 'Position', [cen(1) 10*(max(ZLim)-min(ZLim))+cen(2) cen(3)], ...
        'Style', 'local', 'Color', [1 1 1]/3 );
    light( 'Position', [cen(1) 10*(min(ZLim)-max(ZLim))+cen(2) cen(3)], ...
        'Style', 'local', 'Color', [1 1 1]/3 );
    
end
 
end

%==========================================================================
% CALCULATEBARYCENTRICINTERPOLATIONVALUES
%==========================================================================

function [ lambda1, lambda2, lambda3 ] = ...
    calculateBarycentricInterpolationValues( sizep )
% CALCULATEBARYCENTRICINTERPOLATIONVALUES Calculates the barycentric
% interpolation values of a set of query points in a square of size
% (sizep+1)X(sizep+1) relative to an abstract triangle in the upper region
% of that square

% Define the triangle in the upper part of the square
x1 = sizep; y1 = sizep;
x2 = sizep; y2 = 0;
x3 = 0; y3 = sizep;

% Calculate the barycentric coordinates for each query point on the square
% (These will correspond to the linear indices corresponding to the 2D
% image associated with the triangle)
[x,y] = ndgrid(0:sizep,0:sizep); % NOTE THE USE OF NDGRID
x = x(:); y = y(:);

% detT is MINUS twice the area of the abstract triangle
detT = (x1-x3)*(y2-y3) - (x2-x3)*(y1-y3);
lambda1 = ( (y2-y3).*(x-x3) + (x3-x2).*(y-y3) ) ./ detT;
lambda2 = ( (y3-y1).*(x-x3) + (x1-x3).*(y-y3) ) ./ detT;
lambda3 = 1 - lambda1 - lambda2;

end

