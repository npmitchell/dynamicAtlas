function [xyOut, xyPaths, xyPathsMat] = rk4DynamicVelocity2D(pts, X0v, Y0v, vx, vy, hh, tt, options)
% [xyOut, xyPaths] = rk4DynamicVelocity2D(pts, XY0v, v2d, hh, tt)
% solve Runge-Kutta 4th order in 2d space with fixed velocity field
%
% Parameters
% ----------
% pts : Px2 numeric
%   points at which we begin tracing paths in 2d
% XY0v : MxN numeric
%   evaluation points at which velocities reside
% vx : TxMxN numeric
%   dynamic 2d velocity in x at evaluation points (time, vxgridX, vxgridY)
% vy : TxMxN numeric
%   dynamic 2d velocity in y at evaluation points (time, vygridX, vygridY)
% hh : optional float 0<h<1 (default=0.05)
%   stepsize for finding path
% tt : optional Qx1 float array (default=0:hh:size(v2d, 1))
%   timestamps relative to the indices of v2d at which to compute
%
% Returns
% -------
% xyOut : Px2 float
%   final positions
% xyPaths : length P cell array
%   paths from pts to xyOut
%   
% NPMitchell 2022

if nargin < 8
    options = struct('boundXY', true) ;
    preview = true ;
else
    if ~isfield(options, 'preview')
        preview = true ;
    else
        preview = options.preview ;
    end
    if ~isfield(options, 'boundXY')
        options.boundXY = true ;
    end
    if isfield(options, 'timeAxis')
        timeAxis = options.timeAxis ;
    else
        timeAxis = 1 ;
    end
end

% step size
if nargin < 6
    hh=0.05;                                             
end

% check that v2d is 3D (time, vx, vy)

% parameterization of the path steps in space
if nargin < 7
    tt = 1:hh:size(vx, 1);    
end

% bound the poisitions to lie inside extrapolation region
bound = options.boundXY ;
if bound
    Xmax = max(X0v(:)) ;
    Ymax = max(Y0v(:)) ;
    Xmin = min(X0v(:)) ;
    Ymin = min(Y0v(:)) ;
end
    
% initial condition is xy0

% Handle axes
if timeAxis == 1
    % no additional processing needed
elseif timeAxis == 3
    vx = permute(vx, [3, 1, 2]) ;
    vy = permute(vy, [3, 1, 2]) ;
end

% Preallocation
xyOut = zeros(size(pts)) ;
xyPaths = cell(size(pts, 1), 1) ;
xyPathsMat = zeros(size(pts, 1), length(tt), 2) ;

% calculation loop
for ii=1:(length(tt)-1)         

    disp(['integrating points to time ' num2str(tt(ii)) '/' num2str(max(tt))])
    % for qq = 1:size(pts, 1)
    
    % initial position
    if ii == 1
        x = pts(:, 1) ;
        y = pts(:, 2) ;
        
        % place starting positions in matrix of paths
        if nargout > 1
            for qq = 1:size(pts, 1)    
                xyPathsMat(:, ii, :) = [x,y] ;
            end
        end
    end
    
    % Figure out which time indices to interpolate between for this
    % sub-step ii.
    %
    %          idx0               idx1
    %  ---------o------------------o------------------o---------
    %                ^      ^
    %                |      |
    %              tt(ii)  tt(ii+1)
    %    
    %
    %  TODO: handle this case:
    %
    %          idx0               idx1
    %  ---------o------------------o------------------o---------
    %                          ^      ^
    %                          |      |
    %                        tt(ii)  tt(ii+1)
    %    
    %

    idx0a = floor(tt(ii)) ;
    idx0b = ceil(tt(ii)) ;
    idx1a = floor(tt(ii)+0.5*hh) ;
    idx1b = ceil(tt(ii)+0.5*hh) ;
    idx2a = floor(tt(ii)+hh) ;
    idx2b = ceil(tt(ii)+hh) ;

    % interpolate v2d between time indices idx0 and idx1 using 
    % evaluation points XY0v
    frac1 = 1-(tt(ii)-idx0a) ; % weight for idx0 at the beginning of this sub-step
    halfT = 0.5 * (tt(ii) + tt(ii+1)) ;
    frac2 = 1-(halfT - idx1a) ;
    frac3 = frac2 ;
    frac4 = 1-(tt(ii+1)-idx2a) ; % weight for idx0 by the end of this sub-step

    vx0a = squeeze(vx(idx0a, :, :)) ;
    vx0b = squeeze(vx(idx0b, :, :)) ;
    vx1a = squeeze(vx(idx1a, :, :)) ;
    vx1b = squeeze(vx(idx1b, :, :)) ;
    vx2a = squeeze(vx(idx2a, :, :)) ;
    vx2b = squeeze(vx(idx2b, :, :)) ;
    vy0a = squeeze(vy(idx0a, :, :)) ;
    vy0b = squeeze(vy(idx0b, :, :)) ;
    vy1a = squeeze(vy(idx1a, :, :)) ;
    vy1b = squeeze(vy(idx1b, :, :)) ;
    vy2a = squeeze(vy(idx2a, :, :)) ;
    vy2b = squeeze(vy(idx2b, :, :)) ;
        
    k1x = interp2(X0v, Y0v, ...
        frac1 * vx0a + (1-frac1) * vx0b, x, y, 'makima') ;
    k1y = interp2(X0v, Y0v, ...
        frac1 * vy0a + (1-frac1) * vy0b, x, y, 'makima') ;
    k2x = interp2(X0v, Y0v, ...
        frac2 * vx1a + (1-frac2) * vx1b, ...
        x + 0.5 * hh * k1x, y + 0.5 * hh * k1y, 'makima') ;
    k2y = interp2(X0v, Y0v, ...
        frac2 * vy1a + (1-frac2) * vy1b, ...
        x + 0.5 * hh * k1x, y + 0.5 * hh * k1y, 'makima') ;
    k3x = interp2(X0v, Y0v, ...
        frac3 * vx1a + (1-frac3) * vx1b, ...
        x + 0.5 * hh * k2x, y + 0.5 * hh * k2y, 'makima') ;
    k3y = interp2(X0v, Y0v, ...
        frac3 * vy1a + (1-frac3) * vy1b, ...
        x + 0.5 * hh * k2x, y + 0.5 * hh * k2y, 'makima') ;
    k4x = interp2(X0v, Y0v, ...
        frac4 * vx2a + (1-frac4) * vx2b, ...
        x + k3x * hh, y + k3y * hh, 'makima') ;
    k4y = interp2(X0v, Y0v, ...
        frac4 * vy2a + (1-frac4) * vy2b, ...
        x + k3x * hh, y + k3y * hh, 'makima') ;

    % main equation
    x = x + (1/6)*(k1x + 2*k2x + 2*k3x + k4x)*hh;  
    y = y + (1/6)*(k1y + 2*k2y + 2*k3y + k4y)*hh;  
    
    if bound
        x(x>Xmax) = Xmax ;
        y(y>Ymax) = Ymax ;
        x(x<Xmin) = Xmin ;
        y(y<Ymin) = Ymin ;
        
    end
    
    % final position
    xyOut = [x, y] ;
    
    % Keep track of paths
    if nargout > 1
        for qq = 1:size(pts, 1)    
            xyPathsMat(:, ii+1, :) = xyOut ;
        end
    end
    
    if preview
        clf; quiver(X0v, Y0v, vx0a, vy0a, 0); hold on; plot(x, y, '.')
        xlim([min(X0v(:)), max(X0v(:))])
        ylim([min(Y0v(:)), max(Y0v(:))])
        title(['t = ' num2str(tt(ii))])
        pause(0.01)
    end
end

% Reformat paths to be cell array where each entry is a 2D path
if nargout > 1
    for qq = 1:size(pts, 1)    
        xyPaths{qq} = squeeze(xyPathsMat(qq, :, :)) ;
    end
end


