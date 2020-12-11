function [h1, h2, h3, ax, cax, ax3] = vectorFieldHeatPhaseOnImage(im, ...
    xyfstruct, vx, vy, vscale, options)
%VECTORFIELDHEATPHASEONIMAGE(im, xyfstruct, vx, vy, vscale, options)
%   Plot a vector field (vx,vy) evaluated at grid[xx, yy] on an image im
%
% Parameters
% ----------
% im : PxQ numeric array
%   RGB or grayscale image
% xyfstruct : struct with fields 
%   x : Nx1 or (N*M)x1 or NxM float array (required if no field v)
%       x coordinates of vector field
%   y : Mx1 or (N*M)x1 or NxM float array (required if no field v)
%       y coordinates of vector field
%   v : (N*M)x2 float array (required if no fields x,y)
%       x and y coordinates of vector field
%   f : optional #faces x 3 int array
%       connectivity list of the mesh to color over the image, use only if
%       data is not in 
% vx : N*M x 1 float array
%   velocity in x direction
% vy : N*M x 1 float array
%   velocity in y direction
% vscale : float
%   magnitude associated with maximum color/intensity in velocity image
% qopts : struct with fields
%   outfn : str
%       path to save image if given
%   label : str (default='$|v|$ [$\mu$m / min]')
%       colorbar label. Default is '$|v|$ [$\mu$m / min]' 
%   qsubsample : int (default=10)
%       subsampling factor of the quiver field
%   subsamplingMethod : str ('farthestPoint' or 'random')
%       if nPts > 0 and data is not a structured grid, then subsample the 
%       vector field according to this method. 
%       FarthestPoint is slow but gives approximately equally spaced 
%       vectors in the plane.
%   overlay_quiver : bool (default=true)
%       whether to show the quiverplot overlay
%   qscale : float
%       overall scale of the quivers
%   outfn : str
%       output filename for figure as png 
%   xlim : [minx, maxx]
%       minimum and maximum x values in main axis
%   ylim : [miny, maxy]
%       minimimum and maximum y values in main axis
%
% Returns
% -------
% h1 : handle for imshow
% h2 : handle for imagesc
% h3 : handle for quiverplot
%
% Example usage
% -------------
% % Create example image
% im = peaks ;
% [xxs, yys] = meshgrid(1:size(im, 1), 1:size(im, 2)) ;
% % Create example velocity field on the same field (here xxv=xxs, yyv=yys)
% vx = im .* xxs * 0.02 ;
% vy = im .* yys * 0.01 ;
% sf = sqrt(vx.^2 + vy.^2) ;        % scalar field is velocity magnitude
% options.angle = atan2(vy, vx) ;  % polar field
% % Supply 1d x and y lists
% xyf.x = xxs(1, :) ;
% xyf.y = yys(:, 1)' ;
% vscale = max(abs(sf(:))) ;      % color limit for opacity
% options.visibility = 'on' ;
% options.overlay_quiver = false ;
% vectorFieldHeatPhaseOnImage(rand(size(im)), xyf, vx, vy, vscale, options) ;
%
% See also
% --------
% more limited function: scalarFieldOnImage.m
% less limited but unwieldy: scalarVectorFieldsOnImage.m
%
% NPMitchell 2020

% Default options
labelstr = '$|v|$ [$\mu$m / min]' ;
overlay_quiver = true ;
qsubsample = 5 ;
nPts = 0 ;
subsamplingMethod = 'farthestPoint' ;  % ('farthestPoint' 'random' 'custom')
qscale = 5 ;
quiver_vecfield = [] ;
visibility = 'off' ;

% Unpack options
if nargin > 5
    if isfield(options, 'label')
        labelstr = options.label ;
    end
    if isfield(options, 'overlay_quiver')
        overlay_quiver = options.overlay_quiver ;
    end
    if isfield(options, 'qsubsample')
        qsubsample = options.qsubsample ;
    end
    if isfield(options, 'nPts')
        nPts = options.nPts ;
    end
    if isfield(options, 'qscale') 
        qscale = options.qscale ;
    end
    if isfield(options, 'quiver_vecfield') 
        quiver_vecfield = options.quiver_vecfield ;
    end
    if isfield(options, 'visibility')
        if strcmpi(options.visibility, 'on') || ...
                strcmpi(options.visibility, 'off')
            visibility = lower(options.visibility) ;
        end
    end
    if isfield(options, 'fig')
        figure(options.fig)
    else
        close all
        fig = figure('units', 'normalized', ...
            'outerposition', [0 0 1 1], 'visible', visibility) ;
    end
    if isfield(options, 'subsamplingMethod')
        subsamplingMethod = options.subsamplingMethod ;
    end
else
    options = struct() ;
end

%% Set up the figure
% If grayscale image is passed, convert to RGB
if length(size(im)) < 3
    im = cat(3, im, im, im) ;
end
h1 = imshow(im) ;
hold on;

%% Unpack xystruct. If faces are supplied, then we will use patch
try
    try
        xx = xyfstruct.x ;
        yy = xyfstruct.y ;
    catch
        xx = xyfstruct.v(:, 1) ;
        yy = xyfstruct.v(:, 2) ;
    end
catch
    error('Must supply either fields x,y or field v to xyfstruct')
end
vangle = mod(atan2(vy, -vx), 2* pi) ;
speed = vecnorm([vx(:), vy(:)], 2, 2) ;
if size(speed, 1) == numel(xx) && size(speed, 2) == numel(yy)
    gridded_data = true ;
    ww = length(xx) ;
    hh = length(yy) ;
    speed = reshape(speed, [hh, ww]);

    % Compute angle of the velocity vector
    if ~all(size(vangle) == [hh, ww])
        vangle = reshape(vangle, [hh, ww]) ;
    end

    h2 = imagesc(xx, yy, vangle) ;
    set(h2, 'AlphaData', speed / vscale)
elseif isfield(xyfstruct, 'f') || isfield(xyfstruct, 'faces') 
    gridded_data = false ;
    if isfield(xyfstruct, 'f')
        ff = xyfstruct.f ;
    elseif isfield(xyfstruct, 'faces')
        ff = xyfstruct.faces ;
    end
    if numel(speed) == numel(xx)
        h2 = patch( 'Faces', ff, 'Vertices', [xx(:), yy(:)], ...
            'FaceVertexCData', vangle, 'FaceColor', 'flat', ...
            'EdgeColor', 'none', 'FaceVertexAlphaData', speed / vscale, ...
            'FaceAlpha', 'interp') ;
    elseif length(speed) == size(ff, 1) 
        error('Consider case of speeds defined on faces instead of vertices')
    else
        error('Vector field size does not match coordinates')
    end
elseif any(size(speed) == 1) && ~any(size(xx)==1) && all(size(xx)==size(yy))
    % Try reshaping speed
    speed = reshape(speed, size(xx)) ;
    xx = xx(1,:)  ;
    yy = yy(:,1)' ;
    ww = length(xx) ;
    hh = length(yy) ;
    gridded_data = true ;
    % Compute angle of the velocity vector
    if ~all(size(vangle) == [ww hh])
        vangle = reshape(vangle, [ww hh]) ;
    end
    h2 = imagesc(xx, yy, vangle) ;
    set(h2, 'AlphaData', speed / vscale)
elseif any(size(speed) == 1) && any(size(xx)==1) && any(size(yy)==1)
    % Try reshaping speed to #xx x #yy
    ww = length(xx) ;
    hh = length(yy) ;
    speed = reshape(speed, [ww hh]) ;
    gridded_data = true ;
    % Compute angle of the velocity vector
    if ~all(size(vangle) == [ww hh])
        vangle = reshape(vangle, [ww hh]) ;
    end
    h2 = imagesc(xx, yy, vangle) ;
    set(h2, 'AlphaData', speed / vscale) ;
else
    error('Could not identify data input orientation')
end

ax = gca() ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUIVER 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if overlay_quiver    
    if isempty(quiver_vecfield)
        qvx = vx ;
        qvy = vy ;
    else
        qvx = quiver_vecfield(:, 1) ;
        qvy = quiver_vecfield(:, 2) ;
    end
    if gridded_data
        disp('Taking grid subsampling since data is gridded')
        qvx = reshape(qvx, [hh, ww]) ;
        qvy = reshape(qvy, [hh, ww]) ;
        QX = imresize(qvx, [hh / qsubsample, ww / qsubsample], 'bicubic') ;
        QY = imresize(qvy, [hh / qsubsample, ww / qsubsample], 'bicubic') ;
        xq = 1:qsubsample:ww ;
        yq = 1:qsubsample:hh ;
        [xg, yg] = meshgrid(xx(xq), yy(yq)) ;
    elseif nPts == 0
        disp('Taking linear subsampling since nPts == 0')
        rx = xx(:) ;
        ry = yy(:) ;
        xg = rx(1:qsubsample:end) ;
        yg = ry(1:qsubsample:end) ;
        rvx = qvx(:) ;
        rvy = qvy(:) ;
        QX = rvx(1:qsubsample:end) ;
        QY = rvy(1:qsubsample:end) ;
    else
        if strcmp(subsamplingMethod, 'farthestPoint')
            disp('Subsampling via farthestPoint')
            findpt = [mean(xx(:)), mean(yy(:))] ;
            seedIDx = pointMatch(findpt, [xx(:), yy(:)]) ;
            fpvsOpts = struct() ;
            fpvsOpts.preview = true ;
            sampleIDx = farthestPointVertexSampling(nPts, [xx(:), yy(:)], ...
                ff, seedIDx, fpvsOpts) ;
        elseif strcmp(subsamplingMethod, 'random')
            disp('Subsampling via random sampling')
            sampleIDx = round(random(nPts) * length(xx(:))) ;
        elseif strcmp(subsamplingMethod, 'custom')
            try
                sampleIDx = options.sampleIDx ;
            catch
                error('With custom subsampling method, must supply sampleIDx')
            end
        end
        rx = xx(:) ;
        ry = yy(:) ;
        xg = rx(sampleIDx) ;
        yg = ry(sampleIDx) ;
        rvx = qvx(:) ;
        rvy = qvy(:) ;
        QX = rvx(sampleIDx) ;
        QY = rvy(sampleIDx) ;
    end
    
    h3 = quiver(xg(:), yg(:), qscale * QX(:), qscale * QY(:), 0, ...
        'k', 'LineWidth', 1.2) ;
else
    h3 = [] ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phasemap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap('phasemap') ;
caxis([0, 2*pi]) ;
if isfield(options, 'ylim')
    ylim(options.ylim)  
end
if isfield(options, 'xlim')
    xlim(options.xlim)  
end
if isfield(options, 'title')
    title(options.title, 'Interpreter', 'Latex')
end
set(gca, 'Position', [0 0.11 0.85 0.8]) ;
% Add phasebar
phasebar('location', [0.87, 0.7, 0.1, 0.1]) ;
ax2 = gca() ;
% Add colorbar
cax = axes('Position',[.9 .3 .02 .3]) ;
[~, yyq] = meshgrid(0:4, 0:100) ;
imshow(fliplr(yyq/max(yyq(:)))) ;
axis on
yyaxis right
ylabel(labelstr, 'color', 'k', 'Interpreter', 'Latex') ;
yticks([0 1])
yticklabels({'0', num2str(vscale)})
xticks([])
yyaxis left
yticks([])
cax.YAxis(1).Color = 'k';
cax.YAxis(2).Color = 'k';

% folds
% plot([foldx; foldx], [0, 0, 0; yesz, yesz, yesz], 'k--')
if isfield(options, 'outfn')
    disp(['Saving figure: ', options.outfn])
    saveas(fig, options.outfn) ;   
    close all
end
