function [chisq, chisqn, ssr, ssr_raw, shifts] = chisquareCurves(curv, refcurvsX, ...
    refcurvsY, refvariance, LperiodicX, LperiodicY, smooth_var, optimize_translation, preview)
% chisquareCurves(curv, refcurvs, refvar, take_mean, symmetrize)
% 
% todo: enable refcurvs and re.fvar to be cell instead of 2D array
% 
% Parameters
% ----------
% refcurvsX : Nx1 or NxQ float or int array
%   1d or 2d array of x values for each curv in refcurvs. If 1d, then each
%   reference curve (refcurv) is treated as having the same domain.
% refcurvsY : N x 1 float array
%   The y coordinates of the reference curves. Each row refcurvsY(ii, :) is
%   an array of y coordinates to compare to curv to generate a chisquared
% refvars : N x 1 float array
%   variance of the reference curve at each evaluation site (refcurvsX)
% LperiodicX : length 1 float
%   if first dimension is periodic, give length for modulo
% LperiodicY : length 1 float
%   if second dimension is periodic, give height for modulo
% smooth_var : positive definite or false
%   fractional width along curve for which to smooth the variance using
%   moving median filter
%   
%
% Returns
% -------
% 
%
% NPMitchell 2020

% optional parameters
if nargin < 6
    optimize_translation = false ;
end

% create chisq array
chisq = zeros(size(refcurvsY, 1), 1) ;
if nargout > 1
    chisqn = zeros(size(refcurvsY, 1), 1) ;
    if nargout > 2
        ssr = zeros(size(refcurvsY, 1), 1) ;
        if nargout > 3
            ssr_raw = zeros(size(refcurvsY, 1), 1) ;
        end
    end
end


% compute a chisquare for each reference curve
for ii = 1:size(refcurvsY, 1)
    if all(size(refcurvsX) > 1)
        refcurvX = refcurvsX(ii, :) ;
    else
        refcurvX = refcurvsX ;
    end
    refcurvY = refcurvsY(ii, :) ;
    
    refcurv = [ refcurvX(:), refcurvY(:) ] ;
    refvar = refvariance(ii, :) ;
    
    % mask zeros in variance
    refvar(refvar == 0) = NaN ;
    
    % Optimize the placement of the curve wrt reference
    if optimize_translation
        % disp('chisquareCurves: Optimizing translation')
        % Optimize for DV motion, add to list of leading for this TP
        x0 = [0., 0.] ;
        take_mean = true ;
        symmetrize = true ;
        fun = @(x)ssrCurves(curv + [x(1) x(2)], refcurv, ...
            take_mean, symmetrize) ;
        options = optimset('TolX', 1e-3, 'TolFun', 0.00001);
        shifts = fminsearch(fun, x0, options) ;
        disp(['optimal translation = ', num2str(shifts)])
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if preview
            clf
            subplot(1, 2, 1)
            if smooth_var
                var2plot = movmedian(refvar, smooth_var * length(refvar), 'omitnan') ;
            else
                var2plot = refvar ;
            end
            errorbar(refcurv(:, 1), refcurv(:, 2), [], [], sqrt(var2plot), sqrt(var2plot))
            hold on;
            plot(curv(:, 1), curv(:, 2), '.') ;
            plot(curv(:, 1)+ shifts(1), curv(:, 2)+shifts(2), '.')
            legend({'reference', 'curve', 'optimal curve'})
            hold off;
            axis equal
            subplot(1, 2, 2)
            plot(chisq, '.-')
            if nargout > 1
                yyaxis right
                plot(chisqn, '.-')
            end
            sgtitle(['t= ' num2str(ii) '/' num2str(size(refcurvsY, 1))])
            
            pause(1e-7)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Get raw SSR if requested
        if nargout > 3
            ssr_raw_ii = ssrCurves(curv, refcurv, take_mean, symmetrize) ;
        end
        
        % Now offset the curve for optimization
        curv = curv + [shifts(1), shifts(2)] ;
        
        % Modulo for periodic AP
        if LperiodicX
            curv(:, 1) = mod(curv(:, 1), LperiodicX) ;
        end
        % Modulo for periodic DV
        if LperiodicY
            curv(:, 2) = mod(curv(:, 2), LperiodicY) ;
        end
    else
        % no need to optimize, so ssr will be the same as ssr_raw
        if nargout > 3
            ssr_raw_ii = ssrCurves(curv, refcurv, take_mean, symmetrize) ;
        end
    end    
    
    d2 = pdist2(curv, refcurv, 'squaredeuclidean');
    [dists, idx] = nanmin(d2, [], 2);
    
    % Should we smooth the variance over the whole curve?
    if smooth_var
        % yes, smooth the variance by taking mean variance
        chisq(ii) = nansum(dists' ./ movmedian(refvar(idx), smooth_var * length(idx), 'omitnan'), 2) ;
    else
        % no, use position-dependent variance (may give some nans if
        % variance goes to zero anywhere)
        chisq(ii) = nansum(dists' ./ refvar(idx), 2) ;
    end
    
    
    % Debug 
    % close all 
    % figure ;
    % plot(curv(:, 1), curv(:, 2), '.'); hold on;
    % plot(refcurv(:, 1), refcurv(:, 2), '.') 
    % title('checking that optimize translation worked')
    % waitfor(gcf)
    
    % check result
    % imagesc(d2) 
    % ylabel('curv index')
    % xlabel('refcurv index')
    
    if nargout > 1
        chisqn(ii) = chisq(ii) / sum(~isnan(refvar(idx))) ;
        if nargout > 2
            ssr(ii) = nansum(dists, 1) ;
            
            if nargout > 3
                ssr_raw(ii) = ssr_raw_ii ;
                % we request a raw ssr output. Check that this is the same
                % as ssr just computed if there was no optimization
                if ~optimize_translation
                    assert(all(ssr(ii) == ssr_raw_ii))
                end
            end
        end
    end
end

