function [xyOut, xyPaths] = rk4FixedVelocity2D(pts, X0v, Y0v, vx, vy, hh, tt)
%[xyOut, xyPaths] = rk4FixedVelocity2D(pts, XY0v, v2d, hh, tt)
% solve Runge-Kutta 4th order in 2d space with fixed velocity field
%
% Parameters
% ----------
% pts : Px2 numeric
%   points at which we begin tracing paths in 2d
% XY0v : MxN numeric
%   evaluation points at which velocities reside
% v2d : MxN numeric
%   2d velocity at evaluation points
% hh : float 0<h< 1
%   stepsize for finding path
%
% Returns
% -------
% xyOut : Px2 float
%   final positions
% xyPaths : length P cell array, such that xyPaths{i} is #ttx2 float array
%   paths from pts to xyOut
%   
% NPMitchell 2022

% step size
if nargin < 6
    hh=0.05;                                             
end

% parameterization of the path steps in space
if nargin < 7
    tt = 0:hh:1;    
end

% initial condition is xy0
% interpolate v2d from XY0v
fx = griddedInterpolant(X0v, Y0v,  vx) ;
fy = griddedInterpolant(X0v, Y0v,  vy) ;

% calculation loop
xyOut = zeros(size(pts)) ;
xyPaths = cell(size(pts, 1), 1) ;
for qq = 1:size(pts, 1)
    if mod(qq, 50) == 1 
        disp(['tracing path for pt: ' num2str(qq) '/' num2str(size(pts,1))])
    end
    
    % initial position
    x = pts(qq, 1) ;
    y = pts(qq, 2) ;
    
    % preallocate
    xy = zeros(length(tt), 2) ;
    xy(1, :) = [x, y] ;
    
    for i=1:(length(tt)-1)                              
        % k1: tt(i)
        k1x = fx(x, y);
        k1y = fy(x, y);
        % k2: tt(i)+0.5*hh
        k2x = fx(x + 0.5 * hh * k1x, y + 0.5 * hh * k1y);
        k2y = fy(x + 0.5 * hh * k1x, y + 0.5 * hh * k1y);
        % k3: tt(i)+0.5*hh
        k3x = fx(x + 0.5 * hh * k2x, y + 0.5 * hh * k2y);
        k3y = fy(x + 0.5 * hh * k2x, y + 0.5 * hh * k2y);
        % k4: tt(i)+hh
        k4x = fx((x + k3x * hh), y + k3y * hh);
        k4y = fy((x + k3x * hh), y + k3y * hh);
        % main equation
        x = x + (1/6)*(k1x + 2*k2x + 2*k3x + k4x) * hh;  
        y = y + (1/6)*(k1y + 2*k2y + 2*k3y + k4y) * hh;  
        xy(i+1,:) = [x, y] ;
    end
    
    % final position
    xyOut(qq, :) = xy(length(tt), :) ;
    if nargout > 1
        xyPaths{qq} = xy ;
    end
end


