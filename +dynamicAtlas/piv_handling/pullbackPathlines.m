function [XX, YY] = pullbackPathlines(pivStack, x0, y0, t0, options)
% pullbackPathlines(pivStack, x0, y0, t0, options)
%   Create paths in pullback space (in pixels, XY) by following optical
%   flow of PIV measured on standardized pullback map.
%   Note: I've chosen to use spatial smoothing via Gaussian blur rather
%   than temporal smoothing of velocities for this code.
%
% Parameters
% ----------
% pivStack : struct with two #TP x nX x nY float entries
%   vx: #TP x nX x nY float
%       velocities in x direction
%   vy: #TP x nX x nY float
%       velocities in y direction
%   x : nX x nY numeric
%       piv evaluation coordinates in x
%   y : nX x nY numeric
%       piv evaluation coordinates in y
% x0 : n*m float array 
%   x coordinates in pullback pixels to start pathlines at t0
% y0 : n*m float array 
%   y coordinates in pullback pixels to start pathlines at t0
% t0 : 
%   time at which to begin the pathlines, if member of
%   options.timePoints. Otherwise, treated as an index into the PIV arrays
% options : struct with fields 
%   preview : bool (default=false)
%       view intermediate results
%   timePoints : 1d int array (default= 1:size(pivStack.vx, 1))
%       the timpepoints along which to map pathlines that intersect with
%       (x0, y0) at t0
% 
% Returns
% -------
% XX : #timePoints x size(x0, 1) x size(x0, 2) float array
%   the pullback pixel X coordinates of the pathline spanning timePoints
% YY : #timePoints x size(x0, 1) x size(x0, 2) float array
%   the pullback pixel Y coordinates of the pathline spanning timePoints
%
% NPMitchell 2021

    %% Input checking
    assert(all(size(x0) == size(y0))) 

    %% Default options
    preview = false ;
    debug = false ;

    %% Unpack options
    if isfield(options, 'preview')
        preview = options.preview ;
    end
    if isfield(options, 'debug')
        debug = options.debug ;
    end
    if isfield(options, 'timePoints')
        timePoints = options.timePoints ;
    else
        timePoints = 1:size(pivStack.vx, 1) ;
    end
    % Note: doesn't actually matter whether pullback is a single or double cover.
    % doubleCovered = true ;

    %% Set it up
    first = true ;
    ntps = length(timePoints) ;

    %% First load raw PIV and store (uu, vv) values in a grid
    % Load raw PIV for (uu, vv) in pix/dt
    % Load the positions of the velocity vectors in pixels

    if isfield(options, 'Lx')
        disp('pullbackPathlines(): using supplied Lx')
        Lx = options.Lx ;
    else
        % standard size for pullbacks in atlas
        Lx = 696 ;
    end
    if isfield(options, 'Ly')
        disp('pullbackPathlines(): using supplied Ly')
        Ly = options.Ly ;
    else
        % standard size for pullbacks in atlas
        Ly = 820 ;
    end
    disp('Building pathlines')

    %% Collate velocities into 4d array -- already assumed to be done
    if debug
        % Template example for debugging
        % Debug this function: make fake PIV
        x0 = pivStack.vx{1} ;
        y0 = pivStack.vy{1} ;
        vPIV = zeros(ntps, size(x0, 1), size(x0, 2), 2);
        for ii = 1:ntps
            vPIV(ii, :, :, 1) = 100 * cos(2 * pi * ii / 20) + 5 ;
            vPIV(ii, :, :, 2) = 10 ;
        end
    else
        % Define x and y evaluation points
        xpiv = pivStack.x ;
        ypiv = pivStack.y ;

        % Load in true PIV
        first = true ;
        for ii = 1:ntps
            % tidx = QS.xp.tIdx(timePoints(ii)) ; 
            uu = squeeze(pivStack.vx(ii, :, :)) ;
            vv = squeeze(pivStack.vy(ii, :, :)) ; 

            % Ensure no NaNs in uu and vv
            if any(isnan(uu(:))) || any(isnan(vv(:)))
               disp('inpainting NaNs in uu & vv')
               uu = inpaint_nans(uu) ;
               vv = inpaint_nans(vv) ;
            end

            % Build PIV velocity grid (time, ugrid(1), ugrid(2), x/y)
            if first
                vPIV = zeros(ntps, size(uu, 1), size(uu, 2), 2);
                first = false ;
            end

            vPIV(ii, :, :, 1) = uu ;             % in pix/dt
            vPIV(ii, :, :, 2) = vv ;             % in pix/dt            

        end
    end

    %% Propagate along velocities forward and backward

    % Could use streamline but there is an issue with moving out of the
    % frame. Instead use griddedInterpolant with padded edges, clip at each
    % step along the way to keep in frame, and use periodic BC for vertical
    % direction
    % Preallocate positions for all time
    XX = zeros(length(timePoints), size(x0, 1), size(x0, 2)) ;
    YY = zeros(length(timePoints), size(y0, 1), size(y0, 2)) ;
    % Fill in position at starting time t0
    idx0 = find(timePoints == t0) ;  % QS.xp.tIdx(t0) ;

    % Permit duplicate timePoints for timeaverging at t~0 or t~end, but print 
    % warning when doing so. Do not permit duplicate timePoints in the middle
    % of the array
    if length(idx0) > 1
        disp("WARNING: pullbackPathlines(): multiple timepoints are identical and equal t0")
        for qq = 1:length(idx0)
            ttmp = idx0(qq) ;
            XX(ttmp, :, :) = x0 ;
            YY(ttmp, :, :) = y0 ;
        end
        if ~any(idx0 == 1) && ~any(idx0 == length(timePoints))
            error("WARNING: duplicate points in the middle of timePoint list")
        end
    else
        XX(idx0, :, :) = x0 ;
        YY(idx0, :, :) = y0 ;
    end

    % Now we are done with input handling and timepoint t0.
    if any(diff(timePoints) == 0)
        disp('Error: duplicate timePoints in list -- could allow duplicates at beginning or end...')
    end

    % Propagate forward first: tIdx(t0)+2 onward
    for qq = (max(idx0)+1):length(timePoints)
        disp(['tidx = ' num2str(qq)])
        % 1. Interpolate velocities at time qq-1
        uu = squeeze(vPIV(qq-1, :, :, 1)) ;
        vv = squeeze(vPIV(qq-1, :, :, 2)) ;
        % in transposed coords
        ui = griddedInterpolant(xpiv, ypiv, uu, 'linear', 'nearest') ; 
        vi = griddedInterpolant(xpiv, ypiv, vv, 'linear', 'nearest') ; 

        % 2. Evaluate at XY(qq-1) non-transposed coords
        xx = squeeze(XX(qq-1, :, :)) ;
        yy = squeeze(YY(qq-1, :, :)) ;
        [xclip, yclip] = clipXY(xx, yy, Lx, Ly) ;
        
        % assert(all(abs(xx(:)) > 0))
        dx = reshape(ui(xclip(:), yclip(:)), size(xx)) ;
        dy = reshape(vi(xclip(:), yclip(:)), size(yy)) ;

        % 3. push XY
        Xqq = xx + dx ;
        Yqq = yy + dy ;

        % 4. Clip at x=0,Lx and wrap at y=0,2Ly
        % --> Here we only clip for piv evaluation
        % Clip the coordinates in pullback space
        % [Xqq, Yqq] = clipXY(Xqq, Yqq, Lx, Ly) ;
        
        XX(qq, :, :) = Xqq ;
        YY(qq, :, :) = Yqq ;

        if debug
            subplot(1, 2, 1)
            plot(qq, dx(1, 19, 19), 'b.')
            plot(qq, dy(1, 19, 19), 'r.')
            hold on;
            subplot(1, 2, 2) 
            plot(qq, XX(qq, 19, 19), 'b.')
            plot(qq, YY(qq, 19, 19), 'r.')
            hold on;
        end
    end

    % Propagate backward in time if t0 > timePoints(2)
    if min(idx0) > 1 
        backward_times = fliplr( 1:(min(idx0)-1) ) ;
        for qq = backward_times 
            disp(['tidx = ' num2str(qq)])
            % 1. Interpolate velocities of time qq at their advected 
            %    locations in qq+1.
            %
            %   .---->*     .<----*
            %     qq=.       qq+1=*
            % 
            %    Since the advected x0+u_qq,y0+v_qq are unstructured, we
            %    interpolate current velocity at the advected positions,
            %    ie what velocities got you there, evaluate the velocities 
            %    that pushed to the next positions at the next positions,
            %    and pull back along those velocities.
            %
            uu = squeeze(vPIV(qq, :, :, 1)) ;
            vv = squeeze(vPIV(qq, :, :, 2)) ;

            % Load Lx, Ly from vector (match timepoint)
            [xa, ya] = clipXY(xpiv(:) + uu(:), ypiv(:) + vv(:), Lx, Ly) ;
            ui = scatteredInterpolant(xa, ya, uu(:), 'natural', 'nearest') ;
            vi = scatteredInterpolant(xa, ya, vv(:), 'natural', 'nearest') ;

            % 2. Evaluate at XY(qq+1) non-transposed coords with clipping
            xx = squeeze(XX(qq+1, :, :)) ;
            yy = squeeze(YY(qq+1, :, :)) ;
            [xclip, yclip] = clipXY(xx, yy, Lx, Ly) ;
            

            % 3. Pull XY back
            dx = reshape(ui(xclip(:), yclip(:)), size(xclip)) ;
            dy = reshape(vi(xclip(:), yclip(:)), size(xclip)) ;
            Xqq = xx - dx ;
            Yqq = yy - dy ;

            % 4. Clip at x=0,Lx and wrap at y=0,2Ly
            % [Xqq, Yqq] = clipXY(Xqq, Yqq, Lx, Ly) ;
            
            XX(qq, :, :) = Xqq ;
            YY(qq, :, :) = Yqq ;
        end
    end

    %% Check debugged pathlines
    if debug
        if preview
            for qq = 1:50:size(x0, 1)
                for pp = 1:50:size(x0, 2)
                    plot3(1:length(timePoints), XX(:, qq, pp), YY(:, qq, pp), '-')
                    hold on;
                end
            end
            xlabel('t')
            ylabel('X [pix]')
            zlabel('Y [pix]')
            title('Artificial flow for streamline checking')
        end
        view(2)
        saveas(gcf, fullfile(QS.dir.pivAvg, 'streamline_test_xt.png'))
        view(90, 90)    
        saveas(gcf, fullfile(QS.dir.pivAvg, 'streamline_test_yt.png'))
        view(90, 0)
        saveas(gcf, fullfile(QS.dir.pivAvg, 'streamline_test_xy.png'))
        waitfor(gcf)
    end
end


function [Xqq, Yqq] = clipXY(Xqq, Yqq, Lx, Ly)
    % Note we use minimum values of 1 (in pixels)
    minX = 1 ;
    minY = 1 ;
    
    % % Clip in X
    % Xqq(Xqq > Lx) = Lx ;
    % Xqq(Xqq < minX ) = 1 ;
    % % modulo in Y
    % Yqq(Yqq > Ly) = Yqq(Yqq > Ly) - Ly + minY;
    % Yqq(Yqq < minY) = Yqq(Yqq < minY) + Ly ;

    % Clip in Y
    Xqq(Xqq > Lx) = Xqq(Xqq > Lx) - Lx + minX ;
    Xqq(Xqq < minX) = Xqq(Xqq < minX) + Lx ;
    
    % Modulo in X
    Yqq(Yqq > Ly) = Ly ;
    Yqq(Yqq < minY) = minY ;
    
end