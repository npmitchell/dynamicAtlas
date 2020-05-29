function writeTiffPages(im, fn, options)
% WRTITETIFFPAGES
%   Write a tiff stack to disk page by page
%
% Parameters
% ----------
% im : MxNxP numeric array
%   the numeric array to write to disk as a tiff
% fn : str
%   the filename (full path) to save the tiff
% options : struct, optional
%
%
% NPMitchell 2020

% todo: parse options

for z = 1:size(im, 3)
    % write the first page as overwrite mode, then append
    if z == 1
        imwrite(im(:,:,z), fn, 'tiff', 'Compression','none');
    else
        imwrite(im(:,:,z), fn, 'tiff', 'Compression','none','WriteMode','append');    
    end
end