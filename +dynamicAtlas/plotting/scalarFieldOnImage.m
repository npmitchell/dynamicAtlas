function [h1, h2] = scalarFieldOnImage(im, xx, yy, field, alphaVal, ...
    scale, label, varargin)
%SCALARFIELDONIMAGE(im, xx, yy, field, alphaVal, scale, label)
% Plot a scalar field over an image, colored by magnitude, with const alpha
%
% Parameters
% ----------
% im : 
% xx : N x 1 float array
%   y coordinates of the field evaluation locations
% yy : M x 1 float array
%   y coordinates of the field evaluation locations
% field : NxM float array
%   the scalar field to plot as heatmap
% alphaVal : float
%   the opacity of the heatmap
% scale : float
%   maximum absolute value for the field to be saturated in colormap
% label : str
%   colorbar label, interpreted through Latex by default
% varargin : keyword arguments (optional, default='diverging') 
%   options for the plot, with names
%   'style' : 'diverging' or 'positive' or 'negative'
%   'interpreter' : 'Latex', 'default'/'none'
% 
%
% Returns
% -------
% h1 : handle for imshow
% h2 : handle for imagesc
%
%
% NPMitchell 2020

% Unpack options for style (diverging, positive, negative) and cmap
style = 0 ;     % default is diverging
label_interpreter = 'latex' ; % for colorbar label
if ~isempty(varargin)
    for i = 1:length(varargin)
        if isa(varargin{i},'double') 
            continue;
        end
        if isa(varargin{i},'logical')
            continue;
        end
        if ~isempty(regexp(varargin{i},'^[Ss]tyle','match'))
            stylestr = varargin{i+1} ;
        elseif ~isempty(regexp(stylestr,'^[Ii]nterpreter','match'))
            label_interpreter = varargin{i+1} ;
        end
    end
    if ~isempty(regexp(stylestr,'^[Dd]iverging','match'))
        style = 0 ;
    elseif ~isempty(regexp(stylestr,'^[Pp]ositive','match'))
        style = 1 ;
    elseif ~isempty(regexp(stylestr,'^[Nn]egative','match'))
        style = 2 ;
    end
    
end

% Show the image
h1 = imshow(im) ; hold on;
% Overlay the scalar field
h2 = imagesc(xx, yy, field) ;
alpha(alphaVal) ;
if style == 0
    caxis(gca, [-scale, scale]) ;
    colormap(bwr) ;
elseif style == 1
    caxis(gca, [0, scale]) ;
elseif style == 2
    caxis(gca, [-scale, 0]) ;
end
c = colorbar();
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
c.Label.String = label ;
if ~strcmp(label_interpreter, 'default') && ~strcmp(label_interpreter, 'none')
    c.Label.Interpreter = label_interpreter ;
end
