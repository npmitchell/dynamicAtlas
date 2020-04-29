function [xstar, err] = chisqMinUncertainty(chisq, minN2fit, maxN2fit)


[minchi, idx] = min(chisq) ;

% March away from min until we pass chimin+1

% Find the edge of the chisq above min
pidx = idx ;
thischi = chisq(pidx) ;
while thischi < minchi + 1
    pidx = pidx + 1 ;
    thischi = chisq(pidx) ;
end
pidx = pidx + 1 ;

% Find the edge of the chisq below min
nidx = idx ;
thischi = chisq(nidx) ;
while thischi < minchi + 1
    nidx = nidx - 1 ; 
    thischi = chisq(nidx) ;
end
nidx = nidx - 1;

% Fit the domain of chisquareds to a parabola
if pidx - nidx > minN2fit
    if pidx - nidx  < maxN2fit
        c2fit = chisq(nidx:pidx) ;
    else
        n2shave = floor(0.5 * (maxN2fit - (pidx - nidx))) ;
        id2fit = (nidx + n2shave):(pidx - n2shave) ;
        c2fit = chisq(id2fit) ;
    end
end

%% Fit the data
% fitresult = fit(id2fit, c2fit, 'poly2') ;
% ci = confint(fitresult, 0.95) ;
[pp, SS, mu] = polyfit(id2fit, c2fit, 2) ;
% xtilde = (x - meanx) / stdx
meanx = mu(1) ;
stdx = mu(2) ;
% y = a*((x-m)/s)^2 + b*((x-m)/s) + c ;
aa = pp(1) ;
bb = pp(2) ;
cc = pp(3) ;

% Minimum time is (-b/2a), with minimum y value of -(b^2-4ac) / 4a
minxtilde = -bb / (2* aa) ;
% Convert to xstar via mu transformation
xstar = minxtilde * stdx + meanx ;

% Find where in xstar where y rises by 1





