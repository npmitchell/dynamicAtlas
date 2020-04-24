function [shiftdat, shiftfixed] = shiftImagesX(shiftx, tdat, fixed, fillval)
%SHIFTIMAGESX(shiftx, tdat, fixed) 
%   
%
% Parameters
% ----------
% shiftx : int (or float, dealt as int)
%   the shift applied to tdat
% tdat : N x M float or int image
%   the image to shift relative to fixed
% fixed : N x M float or int image
%   the image to keep fixed but extend in x if necessary
% fillval : float or int
%   value to fill out of bounds pixels
%   default value is to fill with the shift value (to penalize these 
%   pixels in an optimization) 
%
% Returns
% -------
% shiftdat : N+few x M+few image
% shiftfixed : N+few x M+few image
% 
% NPMitchell 2019

% Preallocate the shifted image
if nargin < 4
    shiftdat = zeros(size(fixed, 1), size(fixed, 2) + abs(int16(round(shiftx)))) ;
else
    shiftdat = fillval * ones(size(fixed, 1), size(fixed, 2)+ abs(int16(round(shiftx)))) ;
end

% Create both images
if round(shiftx) > 0
    shiftdat(:, (1+int16(round(shiftx))):end) = tdat ;
    
    % Expanded fixed image
    shiftfixed(size(shiftdat, 1), size(shiftdat,2)) = 0 ;
else
    % Keep tdat on left side of image
    shiftdat(:, 1:(end+int16(round(shiftx)))) = tdat ;
    
    % Expanded fixed image, put fixed data on right side
    newfixed = 0 * shiftdat ;
    newfixed(:, (1-int16(round(shiftx))):end) = fixed ;
    shiftfixed = newfixed ;
end

end

