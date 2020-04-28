function vstat = binnedstats(arr, dx) 
% binnedstats(arr, dim2bin)
% Take statistics of 2nd column of an array based on binning of column 1
% 
% Parameters
% ----------
% arr : 
% dx : float or Qx1 float array
%   spacing between bins (bin width) or edges of bins if array
% 
% Returns
% -------
% vstat
%
% NPMitchell 2020
if nargin < 2
    dx = 1 ;
    edges = (min(arr) - 0.5*dx): dx : (max(arr)+dx*0.5) ;
elseif length(dx) > 1
    edges = dx ;
else
    edges = (min(arr) - 0.5*dx): dx : (max(arr)+dx*0.5) ;
end

% bin the first column
bv = arr(:, 1) ;
[a, EE] = discretize(bv, edges) ;
binc = reshape(0.5 * (EE(1:end-1) + EE(2:end)), [length(EE) - 1, 1]) ;
means = zeros(length(binc), 1) ;
vars = zeros(length(binc), 1) ;
for qq = min(a):max(a)
    % If there are any points in this bin
    if any(a==qq)
        % get the mean of the points in this bin
        means(qq) = mean(arr(a==qq, 2)) ;
        vars(qq) = var(arr(a==qq, 2)) ;
    end
end

vstat = cat(2, binc, means, vars) ;


