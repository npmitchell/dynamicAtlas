function ssdiff = compareCurves(vars, curv1, curv2, tquery)
% Given two curves, interpolate the first in reference to the second and
% find the sum squared differences betweeen the two curves as a function of
% a time offset. This can be used to minimize difference between two time
% offset curves. The second curve should have a shorter domain so that the
% first curve's interpolation covers the seconds' range. 
%
% parameters
% ----------
% vars : a vector containing the variables to be minimized
%     vars1 : float
%       time to shift curv1 relative to curv2's t axis to match the two curves
%     vars2 : float
%       time dilation factor
% curv1 : N x 1 float array
%   first column is independent variable t, second is a value f(t)
% curv1 : N x 1 float array
%   first column is independent variable t, second is a value f(t)
% tquery : M x 1 float array
%   independent variable values at which to interpolate
% 
% Returns 
% -------
% ssdiff : float
%   the sum of the squared differences between curv1(t2) and curv2(t2),
%   where t2 are the timestamps given in curv2
%
% Example usage
% -------------
% 
% 
% NPMitchell 2020

vars1 = vars(1); % time offset
vars2 = vars(2); % time stretch 
% vars3 = vars(3); % y value shift

% Unpack 
t1 = curv1(:, 1) ;
p1 = curv1(:, 2) ;
toff = vars2 * t1 + vars1 ;

t2 = curv2(:, 1) ;
p2 = curv2(:, 2) ;


% Interpolate curv1 at times tquery
curv1intrp = interp1(toff, p1, tquery, 'pchip');
curv2intrp = interp1(t2, p2, tquery, 'pchip');
ssdiff = sum((curv1intrp(tquery) - curv2intrp(tquery)).^2) ;

