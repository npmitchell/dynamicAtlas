function [matchtime, matchtime_unc]  = matchTimeSSR(ssr, nfit)
% matchTimeSSR(ssr, nfit)
% 
% 
% NPMitchell 2020

[~, minID] = min(ssr) ;
wfit = round(nfit*0.5-1) ;
optID = max(1, minID-wfit):min(minID+wfit, length(ssr)) ;
% pentad = ssdsm(optID) ;    
pentad = ssr(optID) ;    

% Fit to parabola
% plot(tstamps, ssds)
% plot(ssds)
% plot(tstamps(optID), ssds(optID))    
% Fit to y = a x^2 + b x + c
[p, S] = polyfit(optID(:), pentad(:), 2) ;
a = p(1) ;
b = p(2) ;
% [y_fit,delta] = polyval(p, trange(pentx) * dt, S);
ci = polyparci(p, S, 0.6827) ;
a_unc = a - ci(1, 1) ;
b_unc = b - ci(1, 2) ;
% Find minimum
if a > 0
    matchtime = -b / (2 * a) ;
    dmdb = -1 / (2 * a) ;
    dmda = 0.5 * b / a^2 ; 
    matchtime_unc = sqrt((dmda * a_unc)^2 + (dmdb * b_unc)^2) ;
else
    % the quadratic fit is so bad that it is upside down
    matchtime = minID ;
    matchtime_unc = NaN ;
end
