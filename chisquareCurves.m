function [chisq, chisqn, ssr] = chisquareCurves(curv, refcurvsX, ...
    refcurvsY, refvars, smooth_var, optimize_translation)
% chisquareCurves(curv, refcurvs, refvar, take_mean, symmetrize)
% 
% todo: enable refcurvs and refvar to be cell instead of 2D array
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
    refvar = refvars(ii, :) ;
    
    % mask zeros in variance
    refvar(refvar == 0) = NaN ;
    
    % Optimize the placement of the curve wrt reference
    if optimize_translation
        % Optimize for DV motion, add to list of leading for this TP
        x0 = [0., 0.] ;
        fun = @(x)ssrCurves(curv + [x(1) x(2)], refcurv, true, true) ;
        shifts = fminsearch(fun, x0) ;
        curv = curv + [shifts(1), shifts(2)] ;
        % Modulo for wDV
        curv(:, 1) = mod(curv(:, 1), 1) ;
    end    
    
    d2 = pdist2(curv, refcurv, 'squaredeuclidean');
    [dists, idx] = nanmin(d2, [], 2);
    if smooth_var
        chisq(ii) = nansum(dists' ./ nanmean(refvar(idx)), 2) ;
    else
        chisq(ii) = nansum(dists' ./ refvar(idx), 2) ;
    end
    
    % check it 
    % imagesc(d2) 
    % ylabel('curv index')
    % xlabel('refcurv index')
    
    if nargout > 1
        chisqn(ii) = chisq(ii) / sum(~isnan(refvar(idx))) ;
        if nargout > 2
            ssr(ii) = nansum(dists, 1) ;
        end
    end
end

