function maxdata = slidingmax(dat, win, rmnans)
% slidingmax(dat, win)
%   sliding one dimensional maximum
%
% Parameters
% ----------
% dat : N x 1 float or int array
%   the data to be "smoothed" via sliding maximum
% win : int
%   window size (full width)
%
% Returns
% -------
% maxdata : N x 1 float or int array
%   the 1d data "smoothed" via sliding max 
% 
% NPMitchell 2020

if nargin < 3
    rmnans = true ;
end

if rmnans
    dat(isnan(dat)) = 0 ;
end

subs = ones(win*length(dat)-(win-1)*win,1);
subs(win+1:win:end) = 2-win;
maxdata = max(reshape(dat(cumsum(subs)),win,[]));

% prepend and postpend maxes smaller than half window range on either end
if mod(win, 2) == 1
    % window size is odd
    beginning = dat(1) * ones(1, floor(win * 0.5)) ; 
    for qq=2:floor(win * 0.5)
        beginning(qq) = max(dat(1:qq)) ;
    end
    ending = dat(end) * ones(1, floor(win * 0.5)) ; 
    for qq=2:floor(win * 0.5)
        ending(qq) = max(dat(end-qq:end)) ;
    end
else
    error("have not coded for even window size. do so here.")
end

maxdata = [beginning maxdata ending ] ;