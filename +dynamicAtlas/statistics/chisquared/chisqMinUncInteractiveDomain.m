function [xstar, err, fit_coefs, nidx, pidx] = ...
    chisqMinUncInteractiveDomain(chisq, minN2fit, maxN2fit, ssr, ...
                                header_preamble, dataStruct)
% chisqMinUncertaintyIntractiveDomain(chisq, minN2fit, maxN2fit)
%
% Parameters
% ----------
% chisq : N x 1 float
%   chi squared values over time
% minN2fit : int
%   minimum number of timepoints to consider for initial fit domain
% maxN2fit : int
%   maximum number of timepoints to consider for initial fit domain
% ssr : optional N x 1 float
%   sum of squared residuals
% header_preamble : optional str
%   string to display in the title
% dataStruct : struct with fields
%   dataframe : QxP numeric array
%       the current frame being timestamped using chisquared
%   refcurves : NxQxP numeric array
%       the reference dynamic data against which the current frame is being
%       timestamped
% 
% Returns
% -------
% xstar : float
%   minimum time of the parabolic fit to chisq
% err : float
%   uncertainty in the minimum of the fit
% fit_coefs : 
%   coefficients of the quadratic fit
%
%
% NPMitchell 2020

if nargin < 6 
    dataframe = [] ;
    refcurves = [] ;
else
    dataframe = dataStruct.dataframe ;
    refcurves = dataStruct.refcurves ;
end
if isempty(dataframe) 
    do_overlay = false ;
else
    do_overlay = true ;
end
if nargin < 5
    header_preamble = '' ;
end
if nargin < 4
    ssr = 0 ;
end

% Grab the minimum
[minchi, idx] = min(chisq) ;

% March away from min until we pass chimin+1
% Find the edge of the chisq above min
pidx = idx ;
thischi = chisq(pidx) ;
while thischi < minchi + 1 && pidx < length(chisq)
    pidx = pidx + 1 ;
    thischi = chisq(pidx) ;
end

% Find the edge of the chisq below min
nidx = idx ;
thischi = chisq(nidx) ;
while thischi < minchi + 1 && nidx > 1
    nidx = nidx - 1 ; 
    thischi = chisq(nidx) ;
end

% Fit the domain of chisquareds to a parabola
if pidx - nidx > minN2fit
    if pidx - nidx > maxN2fit
        while pidx - nidx > maxN2fit
            if pidx - 1 > nidx + 1
                pidx = pidx - 1 ;
            end
            if nidx + 1 < pidx - 1
                nidx = nidx + 1 ;
            end
        end
    end
    id2fit = nidx:pidx ;
else
    while pidx - nidx < minN2fit
        % increase the domain to fit
        if pidx + 1 < length(chisq)
            pidx = pidx + 1 ;
        end
        if nidx - 1 > 0
            nidx = nidx - 1 ;
        end
    end
    id2fit = nidx:pidx ;    
end

% Iteratively update fit domain
fig = figure('visible', 'on') ;
fit_is_fine = false ;
while ~fit_is_fine
    id2fit = nidx:pidx ;  
    c2fit = chisq(id2fit) ;

    %% Fit the data
    % fitresult = fit(id2fit, c2fit, 'poly2') ;
    % ci = confint(fitresult, 0.95) ;
    [coefs, SS] = polyfit(id2fit(:) - idx, c2fit(:), 2) ;
    % % xtilde = (x - meanx) / stdx
    meanx = idx ;
    stdx = 1.0 ;
    % y = a*((x-m)/s)^2 + b*((x-m)/s) + c ;
    aa = coefs(1) ;
    bb = coefs(2) ;
    cc = coefs(3) ;

    fit_coefs.p = coefs ;
    fit_coefs.mu = [meanx stdx] ;
    fit_coefs.S = SS ;
    fit_coefs.ystar = cc - bb^2 / (4*aa) ;

    % Minimum time is (-b/2a), with minimum y value of 1-(b^2-4ac) / 4a
    minxtilde = -bb / (2* aa) ;
    % Convert to xstar via mu transformation
    xstar = minxtilde * stdx + meanx ;

    % Find where in xstar where y rises by 1
    % --> a*xt^2 + b*xt + c = (-(b^2-4ac)) / 4a + 1 = c - b^2/(4a) + 1
    % --> a*xt^2 + b*xt + {c - 1 - (c -b^2/4a) } = 0 
    % --> a*xt^2 + b*xt + {- 1 + [b^2/4a] } = 0 
    % define xt +/- Dxt to be the solutions
    % --> Dxt = sqrt(b^2 - 4*a*cnew) / (2*a) ;
    % --> Dxt = sqrt(b^2 - 4*a*(b^2/4a - 1)) / (2*a) ;
    % --> Dxt = sqrt(b^2 - b^2 + 4a)) / (2*a) ;
    % --> Dxt = sqrt(4a) / (2*a) ;
    % --> Dxt = 1 / sqrt(a) ;
    % Then convert Dxt into Dx = Dxt * stdx ;
    err = stdx / sqrt(aa) ; 
    
    % Is this fit fine?
    % Plot it
    timedense = 1:0.1:length(chisq) ;
    % minimimum of the fit is cstar
    cstar = fit_coefs.ystar ;
    
    % If we have data to overlay, make two subpanels
    if do_overlay
        subplot(1, 2, 1)    
    end
    plot(chisq, '.-')
    hold on;
    yy = polyval(fit_coefs.p, timedense, fit_coefs.S, fit_coefs.mu) ;
    plot(timedense, yy, '--')
    plot(nidx:pidx, chisq(nidx:pidx), 'o')
    if length(ssr) == length(chisq) 
        plot(1:length(chisq), ssr / ssr(idx), 'x')
    end
    errorbar(xstar, cstar, err, 'horizontal')
    ylim([min(cstar * 1.1, 0), max(30, cstar * 2)])
    %xlims = xlim() ;
    xlim([0 length(chisq)])
    ylabel('\chi^2 / N')
    xlabel('time')
    legend({'$\chi^2/N$', 'fit', 'data for fit',...
        '$\Sigma(x-\bar{x})^2$ (normalized)'}, ...
        'Location', 'best', 'Interpreter', 'Latex')
    suptitle([header_preamble 'Enter=OK, up=wider, down=narrower, right/left=shift'])
    
    % Add overlay data to other subplot if supplied
    if do_overlay
        subplot(1, 2, 2)
        %imshow(mat2gray(dataframe))
        dataforfig = double(dataframe);
        imshow(dataforfig,[0,prctile(dataforfig,99,'all')])
        hold on;
        indx = min(max(1, round(xstar)), size(refcurves, 2)) ;
        xcurv = squeeze(refcurves(:, indx, 1)) ;
        ycurv = squeeze(refcurves(:, indx, 2)) ;
        plot(xcurv * size(dataframe, 2), ycurv * size(dataframe, 1), '.')
    end
    
    % Update domain
    move_on = false ;
    while ~move_on
        button = waitforbuttonpress() ;
        if button && strcmp(get(gcf, 'CurrentKey'), 'return')
            move_on = true ;
            fit_is_fine = true ;
        elseif button && strcmp(get(gcf, 'CurrentKey'), 'downarrow')
            move_on = true;
            nidx0 = nidx ;
            nidx = min(pidx - 2, nidx + 1) ;
            pidx = max(nidx0 + 2, pidx - 1) ;
            clf
        elseif button && strcmp(get(gcf, 'CurrentKey'), 'uparrow')
            move_on = true;
            nidx = max(1, nidx - 1) ;
            pidx = min(length(chisq), pidx + 1) ;
            clf
        elseif button && strcmp(get(gcf, 'CurrentKey'), 'rightarrow')
            move_on = true;
            nidx = min(pidx - 2, nidx + 1) ;
            pidx = min(length(chisq), pidx + 1) ;
            clf
        elseif button && strcmp(get(gcf, 'CurrentKey'), 'leftarrow')
            move_on = true;
            nidx = max(1, nidx - 1) ;
            pidx = max(nidx + 2, pidx - 1) ;
            clf
        end
        if ~move_on
            disp('Bad button press, try again: Enter, up, or down arrow')
        end
    end  
end


