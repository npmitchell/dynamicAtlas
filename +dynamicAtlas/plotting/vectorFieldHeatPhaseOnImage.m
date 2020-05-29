function [h1, h2, h3, ax, cax, ax3] = vectorFieldHeatPhaseOnImage(im, xx, yy, vx, vy, vscale, ...
    options)
%VECTORFIELDHEATPHASEONIMAGE(im, xx, yy, vx, vy, vscale, options)
%   Plot a vector field (vx,vy) evaluated at grid[xx, yy] on an image im
%
% xx : N x 1 float array
%   x values of PIV grid evaluation points
% yy : M x 1 float array
%   y values of PIV grid evaluation points
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
% NPMitchell 2020

% Default options
labelstr = '$|v|$ [$\mu$m / min]' ;
overlay_quiver = true ;
qsubsample = 10 ;
qscale = 10 ;

% Unpack options
if isfield(options, 'label')
    labelstr = options.label ;
end
if isfield(options, 'overlay_quiver')
    overlay_quiver = options.overlay_quiver ;
end
if isfield(options, 'qsubsample')
    qsubsample = options.qsubsample ;
end
if isfield(options, 'qscale') 
    qscale = options.qscale ;
end


% 
% vangle = reshape(mod(atan2(vy, -vx), 2* pi), gridsz) ;
% speed = reshape(vecnorm([v2dsm_ii(:, 1), v2dsm_ii(:, 2)], 2, 2), gridsz);
ww = length(xx) ;
hh = length(yy) ;
vangle = mod(atan2(vy, -vx), 2* pi) ;
speed = reshape(vecnorm([vx(:), vy(:)], 2, 2), [hh, ww]);

% Compute angle of the velocity vector
if ~all(size(vangle) == [hh, ww])
    vangle = reshape(vangle, [hh, ww]) ;
end

% Set up the figure
close all
fig = figure('units', 'normalized', ...
    'outerposition', [0 0 1 1], 'visible', 'off') ;
h1 = imshow(im) ;
hold on;
h2 = imagesc(xx, yy, vangle) ;
set(h2, 'AlphaData', speed / vscale)
ax = gca() ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUIVER 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if overlay_quiver    
    vx = reshape(vx, [hh, ww]) ;
    vy = reshape(vy, [hh, ww]) ;
    QX = imresize(vx, [hh / qsubsample, ww / qsubsample], 'bicubic') ;
    QY = imresize(vy, [hh / qsubsample, ww / qsubsample], 'bicubic') ;
    xq = 1:qsubsample:ww ;
    yq = 1:qsubsample:hh ;
    [xg, yg] = meshgrid(xx(xq), yy(yq)) ;

    h3 = quiver(xg(:), yg(:), qscale * QX(:), qscale * QY(:), 0, 'k', 'LineWidth', 1.2) ;
else
    h3 = [] ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phasemap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap phasemap
caxis([0, 2*pi])
if isfield(options, 'ylim')
    ylim(ylim)  % [size(im, 2) * 0.25, size(im, 2) * 0.75]
end
if isfield(options, 'xlim')
    ylim(xlim)  
end
set(gca, 'Position', [0 0.11 0.85 0.8]) ;
% Add phasebar
phasebar('location', [0.87, 0.7, 0.1, 0.1]) ;
ax2 = gca() ;
% Add colorbar
cax = axes('Position',[.9 .3 .02 .3]) ;
[~, yyq] = meshgrid(0:4, 0:100) ;
imshow(fliplr(yyq/max(yyq(:))))
axis on
yyaxis right
ylabel(labelstr, 'color', 'k', ...
    'Interpreter', 'Latex')
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
    saveas(fig, options.outfn) ;   
    close all
end
