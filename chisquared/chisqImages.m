function [chisq] = chisqImages(im1, im2, sigma)
%CHISQIMAGES(im1, im2, sigma) Find ChiSquared for difference between images
%   Find chisq = sum{ 1/sigma_i [y_i - y(x_i)]^2 }
%
% Parameters
% ----------
% im1 : N x M matrix
%   the first image
% im2 : N x M matrix
%   the second image
% sigma : N x M matrix or float
%   The standard deviation at each location, or if uniform, the stdev
%
% Returns 
% -------
% chisq : float
%   The chisquared value for the difference in matrices 
%
% NPMitchell 2019 

seq = (im1 - im2).^2 ./ (sigma.^2) ;
chisq = sum(seq(:)) ;

end

