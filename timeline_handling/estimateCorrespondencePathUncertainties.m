function [uncertainties, uncs, movmeanWindowForUnc] = ...
    estimateCorrespondencePathUncertainties(tpath, tpath0, exptIDi, exptIDj, movmeanWindowForUnc)
%[uncertainties, uncs, movmeanWindowForUnc] = ...
%     estimateCorrespondencePathUncertainties(tpath, tpath0, movmeanWindowForUnc)
%
% Parameters
% ----------
% tpath : N x 2
%   each row: timestamp of dataset i, corresponding timestamp of dataset j
% tpath0 : Nx2 float
%   initial guess for time-time-correpondences
%   each row: timestamp of dataset i, guessed correspondence in dataset j
% exptIDi : string
%   string identifier for experiment i
% exptIDj : string
%   string identifier for experiment j
% movmeanWindowForUnc : int (default=10)
%   window size for moving mean
%
% Returns
% -------
% uncertainties : Nx1 float
%   smoothed, median filtered difference between tpath and tpath0
% uncs : Nx1 float
%   median filtered difference between tpath and tpath0
% movmeanWindowForUnc : int
%   window size for moving mean
% 
% NPMitchell 2022

if nargin < 5
    movmeanWindowForUnc = 10 ;
end
uncOk = false ;
while ~uncOk
    guesses = interp1(tpath0(:, 1), tpath0(:, 2), tpath(:, 1)) ;
    % remove nans using median filter
    uncs = abs(tpath(:, 2) - guesses) ;
    uncs = medfilt1m(uncs, 1) ;
    uncertainties = movmean(uncs, movmeanWindowForUnc) ;

    clf
    % Note: I have debugged this extensively. The
    % labels are now correct. 
    % assert(size(XYs{ii}, 1) == size(cij, 1))
    % assert(size(XYs{jj}, 1) == size(cij, 2))
    % imagesc(1:size(XYs{jj}, 1), 1:size(XYs{ii}, 1), cij)
    imagesc(cij)
    title(['c(i,j): ' exptIDi ' to ' exptIDj])
    ylabel([exptIDi ' timeline i [min]'])
    xlabel([exptIDj ' timeline j [min]'])
    cb = colorbar ;
    ylabel(cb, corr_method)
    caxis([-1,1])
    colormap(bwr)
    hold on;
    errorbar(tpath(:, 2), tpath(:, 1), ...
        [], [],uncertainties, uncertainties,'o-', ...
        'color', yellow)

    axis equal
    axis tight
    set(gca,'fontsize', 12);
    msg = 'Uncertainty smoothing window ok? [Y/n]' ;
    title(msg)
    reply = input(msg, 's') ;
    if contains(lower(reply), 'n')
        % Update movmeanWindowForUncertainties
        movmeanWindowForUnc = input('Enter new value for Uncertainty window smoothing = ') ;
    else
        uncOk = true ;
    end
end