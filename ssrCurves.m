function ssr = ssrCurves(curv, refcurv, take_mean, symmetrize)
% ssrCurves(curv, refcurv, take_mean, symmetrize)
% Find sum of squared residual
%
% Parameters
% ----------
% curv : N x 2 float array
%   the curv to match refcurv
% refcurv : M x 2 float array
%   the second curv, reference
% take_mean : bool (default=false)
%   use mean distance between the two (discrete) curves instead of sum
% symmetrize: bool (default=false)
%   return sqrt(ssr_12)* sqrt(ssr_21) so that curve ordering is irrelevant, 
%   or sum them together if take_mean==false. 
%
% Returns
% -------
% ssr : float
%
% NPMitchell 2020 

if nargin < 3 
    take_mean = false ;
end

d1 = pdist2(curv, refcurv);
if take_mean 
    ssr = nanmean(nanmin(d1, [], 2));
    if symmetrize
        d2 = pdist2(refcurv, curv);
        mindist2 = nanmean(nanmin(d2, [], 2));
        ssr = sqrt(ssr) * sqrt(mindist2) ;
    end
else
    ssr = nansum(nanmin(d1, [], 2));
    if symmetrize
        d2 = pdist2(refcurv, curv);
        mindist2 = nanmin(d2, [], 2);
        ssr = ssr + nansum(mindist2);
    end
end

% returns the sum of squared residuals or mean of squared residuals
    