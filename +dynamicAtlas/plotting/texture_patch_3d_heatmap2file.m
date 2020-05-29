function texture_patch_3d_heatmap2file(FF, VV, TF, TV, IV, tMapOpts, figOpts)
% TEXTURE_PATCH_3D_HEATMAP2FILE(FF, VV, TF, TV, IV, tMapOpts, figOpts)
 % unfinished?
 % NPMitchell 2020

% Default figOpts
labelstr = '' ;
interpreter = 'latex' ;
qscale = 10 ;
sscale = 0 ;
alphaVal = 0.8 ;
style = 'diverging' ;
qsubsample = 1 ;
xlabelstr = '' ;
ylabelstr = '' ;
zlabelstr = '' ;
titlestr = '' ;
% figure parameters
figWidth = 16 ;  % cm
figHeight = 10 ; % cm
lw = 1.2 ;       % linewidth
axposition = [0 0.11 0.85 0.8] ;
axisOff = true ;

% Unpack figOpts
if isfield(figOpts, 'style') 
    style = figOpts.style ;
end
if isfield(figOpts, 'sscale') 
    sscale = figOpts.sscale ;
end
if isfield(figOpts, 'label')
    labelstr = figOpts.label ;
end
if isfield(figOpts, 'title')
    titlestr = figOpts.title ;
end
if isfield(figOpts, 'qscale') 
    qscale = figOpts.qscale ;
end
if isfield(figOpts, 'qsubsample')
    qsubsample = figOpts.qsubsample ;
end
if isfield(figOpts, 'alpha') 
    alphaVal = figOpts.alpha ;
end
if isfield(figOpts, 'figWidth') 
    figWidth = figOpts.figWidth ;
end
if isfield(figOpts, 'figHeight') 
    figHeight = figOpts.figHeight ;
end
if isfield(figOpts, 'linewidth') 
    lw = figOpts.linewidth ;
end
if isfield(figOpts, 'axPosition') 
    axposition = figOpts.axPosition ;
end
if isfield(figOpts, 'axisOff') 
    axisOff = figOpts.axisOff ;
end

% PLOT THE TEXTURE PATCH
texture_patch_3d(FF, VV, TF, TV, IV, tMapOpts) ;

% Add title (optional)
if isfield(figOpts, 'title')
    title(figOpts.title, 'Interpreter', 'latex')
end

% Add the colorbar in the style set in figOpts struct
if strcmp(style, 'phase')
    disp('setting scalar field to phase style')
    %%%%%%%%%%%%%%%%%%%
    % Phasemap
    %%%%%%%%%%%%%%%%%%%
    colormap phasemap
    caxis([0, 2*pi])
    set(gca, 'Position', axposition) ;
    % Add phasebar
    phasebar('location', [0.87, 0.7, 0.1, 0.1]) ;
    % Add colorbar
    cax = axes('Position',[.9 .3 .02 .3]) ;
    [~, yyq] = meshgrid(0:4, 0:100) ;
    imshow(fliplr(yyq/max(yyq(:))))
    axis on
    yyaxis right
    ylabel(labelstr, 'color', 'k', ...
        'Interpreter', interpreter)
    yticks([0 1])
    yticklabels({'0', num2str(sscale)})
    xticks([])
    yyaxis left
    yticks([])
    cax.YAxis(1).Color = 'k';
    cax.YAxis(2).Color = 'k';
elseif strcmp(style, 'diverging')
    disp('setting scalar field to diverging style')
    colormap bwr
    if isfield(figOpts, 'ylim')
        ylim(figOpts.ylim)
    end
    set(gca, 'Position', axposition) ;
    if axisOff
        axis off
    end
    
    % Set color axis limits
    if sscale > 0
        caxis([-sscale, sscale])
    end
    
    % Add colorbar
    c = colorbar('Position',[.9 .333 .02 .333]) ;
    % ylabel(cax, labelstr, 'color', 'k', ...
    %     'Interpreter', interpreter)
    
    % Make colorbar share the alpha of the image
    % Manually flush the event queue and force MATLAB to render the colorbar
    % necessary on some versions
    drawnow
    % Get the color data of the object that correponds to the colorbar
    cdata = c.Face.Texture.CData;
    % Change the 4th channel (alpha channel) to 10% of it's initial value (255)
    cdata(end,:) = uint8(alphaVal * cdata(end,:));
    % Ensure that the display respects the alpha channel
    c.Face.Texture.ColorType = 'truecoloralpha';
    % Update the color data with the new transparency information
    c.Face.Texture.CData = cdata;
    c.Label.Interpreter = interpreter ;
    c.Label.String = labelstr ;

else
    error('have not coded for this style yet')
end


% Save the image if outfn is supplied
if isfield(figOpts, 'outfn')
    disp(['scalarVectorFieldsOnImage: saving ' figOpts.outfn])
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 figWidth figHeight]);  
    saveas(gcf, figOpts.outfn) ;   
    close all
end
