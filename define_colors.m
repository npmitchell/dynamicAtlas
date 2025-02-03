function [colors, names] = define_colors(varargin)
%DEFINE_COLORS Define a set of pleasant colors
%   Define a set of colors to use

if length(varargin) > 1
    ncolors = varargin{1} ;
    disp(['Defining ' num2str(ncolors) ' colors...'])
else
    ncolors = 7 ;
end

blue    = [0.0000, 0.4470, 0.7410] ; % 1
red     = [0.8500, 0.3250, 0.0980] ; % 2
yellow  = [0.9290, 0.6940, 0.1250] ; % 3
purple  = [0.4940, 0.1840, 0.5560] ; % 4 
green   = [0.4660, 0.6740, 0.1880] ; % 5 
sky     = [0.3010, 0.7450, 0.9330] ; % 6 
maroon  = [0.6350, 0.0780, 0.1840] ; % 7

colors = [blue; red; yellow; purple; green; sky; maroon] ;

if ncolors < 8
    colors = colors(1:ncolors, :) ;
else
    error('Need to define more colors for this')
end

if nargout > 1
    names = {'blue', 'red', 'yellow', 'purple', 'green', 'sky', 'maroon'};
    names = names(1:ncolors) ;
end

end

