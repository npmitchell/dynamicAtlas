function [xstar, err, fit_coefs] = chisqMinUncertainty(chisq, minN2fit, maxN2fit)


[minchi, idx] = min(chisq) ;

% March away from min until we pass chimin+1

% Find the edge of the chisq above min
pidx = idx ;
thischi = chisq(pidx) ;
while thischi < minchi + 1
    pidx = pidx + 1 ;
    thischi = chisq(pidx) ;
end

% Find the edge of the chisq below min
nidx = idx ;
thischi = chisq(nidx) ;
while thischi < minchi + 1
    nidx = nidx - 1 ; 
    thischi = chisq(nidx) ;
end

% Fit the domain of chisquareds to a parabola
if pidx - nidx > minN2fit
    if pidx - nidx  < maxN2fit
        id2fit = nidx:pidx ;
    else
        n2shave = floor(0.5 * (maxN2fit - (pidx - nidx))) ;
        id2fit = (nidx + n2shave):(pidx - n2shave) ;
    end
else
    while pidx - nidx < minN2fit
        % increase the domain to fit
        pidx = pidx + 1 ;
        nidx = nidx - 1 ;
    end
    id2fit = nidx:pidx ;    
end
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

if nargout > 2
    fit_coefs.p = coefs ;
    fit_coefs.mu = [meanx stdx] ;
    fit_coefs.S = SS ;
    fit_coefs.ystar = cc - bb^2 / (4*aa) ;
end

% Minimum time is (-b/2a), with minimum y value of 1-(b^2-4ac) / 4a
minxtilde = -bb / (2* aa) ;
% Convert to xstar via mu transformation
xstar = minxtilde * stdx + meanx ;

% Find where in xstar where y rises by 1
% --> a*xt^2 + b*xt + c = (1-(b^2-4ac)) / 4a + 1 
% --> a*xt^2 + b*xt + {c - 1 - [1-(b^2-4ac)] / 4a } = 0 
% define xt +/- Dxt to be the solutions
% --> Dxt = sqrt(b^2 - 4*a* {c - 1 - [1-(b^2-4ac)] / 4a }) / (2*a) ;
% Then convert Dxt into Dx = Dxt * stdx ;
cterm = cc - 1 - ((1-(bb^2 -4*aa*cc)) / (4*aa)) ;
err = stdx * (sqrt(bb^2 - 4*aa* cterm ) / (2*aa)) ; 






