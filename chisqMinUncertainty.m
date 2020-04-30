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






