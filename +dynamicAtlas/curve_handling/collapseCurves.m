function [curv_transf, transf, ssr]  = collapseCurves(curv, refcurv, tquery)
% [curv1, curv2, transf]  = collapseCurves(curv1, curv2)
% Optimize for translation of a curve against curv2
%
% Parameters
% ----------
% curv : Mx2 numeric
%   Curve to be shifted as closely as possible to refcurv
% refcurv : Nx2 numeric
%   curve against which the first curv is collapsed
% tquery : Q x 1 numeric
%   timepoints in fixed (reference) parameterization (t) at which the
%   interpolated, shifted curve is to be evaluated against the reference 
%   curve.
%
% Returns
% -------
% curv_transf
% transf
% ssr

x0 = [0., 0.] ;
fun = @(x)compareCurves([x(1) x(2)], curv, refcurv, tquery) ;
options = optimset('TolX', 1e-3, 'TolFun', 0.001);
transf = fminsearch(fun, x0, options) ;

% Transform the curve (scaling + simple translation)
t1 = curv(:, 1) ;
toff = transf(1) * t1 + transf(2) ;
curv_transf = [toff, curv(:, 2)] ;

% Get raw SSR if requested
if nargout > 3
    curv1intrp = interp1(toff, p1, tquery, 'pchip', 'nearest');
    curv2intrp = interp1(t2, p2, tquery, 'pchip', 'nearest');
    ssr = sum((curv1intrp(tquery) - curv2intrp(tquery)).^2) ;
end
