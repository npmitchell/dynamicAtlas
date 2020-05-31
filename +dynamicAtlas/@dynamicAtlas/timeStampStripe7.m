function timeStampStripe7(da, genotype, label, Options)
% TIMESTAMPSTRIPE7(da, genotype, label, Options)
%   Find appropriate timestamp for data stained by 'label' for all fixed
%   pullbacks compared to master timeline using leading edge of stripe7
%   
% Parameters
% ----------
% da : dynamicAtlas class instance
%   the dynamicAtlas for which we find timestamps
% genotype : str
% stain : str
% Options : struct with optional fields
%   save_fancy : optional bool, default=true
%   overwrite : optional bool, default=false
%       overwrite previous results
%   hands_on : optional bool, default=true
%       use interactive domain selection for chisquare fitting
%   cdf_min : optional float, default=0.01
%       intensity cumulative distribution function minimum 
%       cutoff
%   cdf_max : optional float, default=0.999 
%       intensity cumulative distribution function maximum 
%       cutoff
%   sigma : optional float, default=20 
%       smoothing used in the stripe ID in iLastik
%
%
% Returns
% -------
%
% Outputs
% -------
% dynamicAtlas.path/genotype/label/embryoID/timematch_curve7_chisq.mat
% dynamicAtlas.path/genotype/label/embryoID/timematch_curve7_chisq.txt
%
% NPMitchell 2020

% Unpack Options for what to perform
save_fancy = true; 
overwrite = false ;
hands_on = true ;
cdf_min = 0.01 ;
cdf_max = 0.999 ;
% sigma = smoothing used in the stripe ID in iLastik
sigma = 20 ;

% Unpack Options
if nargin > 3
    if isfield(Options, 'preview')
        save_fancy = Options.save_fancy ;
    end
    if isfield(Options, 'overwrite')
        overwrite = Options.overwrite ;
    end
    if isfield(Options, 'hands_on')
        hands_on = Options.hands_on ;
    end
    if isfield(Options, 'cdf_min')
        cdf_min = Options.cdf_min ;
    end
    if isfield(Options, 'cdf_max')
        cdf_max = Options.cdf_max ;
    end
    if isfield(Options, 'sigma')
        sigma = Options.sigma ;
    end
end

% Unpack more options
prepend = da.lookup(genotype).prepend ;
exten = da.lookup(genotype).exten ;
if strcmp(label, 'Runt')
    timerfn = 'timematch_RuntNanobody_stripe7.mat' ;
elseif strcmp(label, 'Eve')
    timerfn = 'timematch_EveYFP_stripe7.mat' ;
else
    error('Label is not Runt or Eve. What stripe7 is this? Code for it here.')
end
step = 1 ;
sigmastep = sprintf('_sigma%03d_step%03d', sigma, step) ;
dt = 1 ;

% Build reference Directory for where the timeline is stored (with the
% timeline's stripe 7 curves)
refDir = dir(fullfile(da.path, 'timing', genotype, label, [da.timeLineMethod 'corr_ss*'])) ;
refDir = fullfile(refDir(1).folder, refDir(1).name) ;


%% Build the lookuptable for just non-dynamic data (assumes dynamic data) 
lum = da.findStaticGenotypeLabel(genotype, label) ;
% note: similar to lum = da.lookup(genotype).map(label) ;

%% Plotting
[colors, names] = define_colors ;
colors = colors ./ vecnorm(colors, 2, 2) ;
yellow = colors(3, :) ;
green = colors(5, :) ;
sky = colors(6, :) ;
stripecolor = green ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assumes we have already run EnsembleStripe to make smoothed images 

%% Process probabilities in iLastik -- train on stripe7 versus not a stripe

%% Identify stripes if not done so already
stripefn = 'Runt_stripe7curve.mat' ;
probfn = [sigmastep(2:end), filesep, 'smooth', filesep, label, '_smooth', ...
    sigmastep, '_bin_Probabilities.h5'] ;
for kk = 1:length(lum.folders)
    
    filename = lum.names{kk} ;
    embryoDir = lum.folders{kk} ;
    embryoID = lum.embryoIDs{kk} ;
    
    stripefnkk = fullfile(embryoDir, stripefn) ;
    if ~exist(stripefnkk, 'file') || overwrite
        disp(['Building ' stripefnkk])

        pfn = fullfile(embryoDir, probfn) ;
        disp(['Opening ' pfn])
        dat = h5read(pfn, '/exported_data') ;
        midx = round(0.5 * size(dat, 2)) ;
        midy = round(0.5 * size(dat, 3)) ;
        dcrop = squeeze(dat(2, midx:end, :)) ;
        addx = midx ;
        addy = 0 ;

        % Extract the leading curve
        maskfn = fullfile(embryoDir, 'Runt_mask.tif');  
        
        thres = 0.5 ;
        minmaxsz = [5e3, 1e8] ;
        curv = extractStripeEdges(dat, maskfn, thres, minmaxsz, embryoID) ;
        % depth = max(curv) - min(curv) ;
        stripe7curve = curv ; 
        stripe7curve_frac = curv ;
        stripe7curve_frac(:, 1) = double(stripe7curve_frac(:, 1)) / double(size(dat, 2)) ;
        stripe7curve_frac(:, 2) = double(stripe7curve_frac(:, 2)) / double(size(dat, 3)) ;
        save(stripefnkk, 'stripe7curve', 'stripe7curve_frac')
        
        % resize stripe7 curv
        curv = stripe7curve_frac ;
    else
        disp(['Loading stripe7curve from ' stripefnkk])
        load(stripefnkk, 'stripe7curve_frac') ;
        curv = stripe7curve_frac ;
    end

    % Save fancy image
    fancyImFn = fullfile(embryoDir, [label 'stripe7_rgb_' embryoID '.png']) ;
    if save_fancy && (~exist(fancyImFn, 'file') || overwrite ) 
        disp(['Building ' fancyImFn])
        close all
        % Load image and all other images of this embryo
        estruct = da.findEmbryo(embryoID) ;

        nlabels = length(estruct.names) ;
        % Make RGB image
        alllabels = '' ;
        for qq = 1:nlabels
            % load this channel/label
            imq = double(imread(fullfile(estruct.folders{qq}, estruct.names{qq}))) ;
            if qq == 1
                combined = zeros(size(imq, 1), size(imq, 2), 3) ;
            end

            % append to all labels for fancy image saving
            alllabels = [alllabels estruct.labels{qq} '_'] ;

            % Adjust intensity
            [f,x] = ecdf(imq(:));
            f1 = find(f>cdf_min, 1, 'first');
            f2 = find(f<cdf_max, 1, 'last');
            lim = [x(f1) x(f2)];
            imq = mat2gray(imq, double(lim));
            combined(:, :, 1) = combined(:, :, 1) + imq * colors(qq, 1) ;
            combined(:, :, 2) = combined(:, :, 2) + imq * colors(qq, 2) ;
            combined(:, :, 3) = combined(:, :, 3) + imq * colors(qq, 3) ;
        end
        % correct for oversaturation
        combined(combined > 1) = 1.0 ;
        combined = combined / max(combined(:)) ;
        wDV = size(combined, 1) ;
        wAP = size(combined, 2) ;
        
        % Now make the figure
        fig = figure('visible', 'off') ;
        imshow(combined)
        hold on;
        imwrite(combined, fullfile(embryoDir, [alllabels 'rgb_' embryoID '.png']))
        plot(curv(:, 1)*wAP, curv(:, 2)*wDV, 'o-', 'color', stripecolor)
        % plot(curv' + midx, (1:length(curv)) + midy - width, 'o-', 'color', yellow)
        title([embryoID ': Identification of stripe 7'])
        saveas(fig, fancyImFn)
        close all
    end   
end

%% Now match to time domain Runt Nanobody or EveYFP data
close all
optimize_trans = true ;

% Load reference curves
tmp = load(fullfile(refDir, 'curve7stats_collapsed_filtered.mat')) ;
LEX = tmp.LXs ;
LES = tmp.LSs ;
LEY = zeros(size(LEX)) ;
for qq=1:size(LEX, 2)
    LEY(:, qq) = (1:size(LEX, 1)) / size(LEX, 1) ;
end

% Go through each embryo and find time
for kk = 1:length(lum.folders)
    embryoDir = lum.folders{kk} ;
    embryoID = lum.embryoIDs{kk} ;
    eDirs4ID = da.findEmbryo(embryoID) ;

    % assert that the current expt is dynamic
    assert(lum.nTimePoints(kk) == 1)
    
    % Which adjusted curve
    load(fullfile(embryoDir, stripefn), 'stripe7curve_frac') ;
    curv = stripe7curve_frac ;
    
    % Match to dynamic data curv
    ssds = [] ; 
    tstamps = [] ;
    options = optimset('MaxIter', 25, 'TolX', 1e-2) ; 
    % Other options: ('FunctionTolerance', 1e-7, 'Display','iter','PlotFcns',@optimplotfval);
    
    % Consider each timepoint, compare to this curve to reference curves
    msg = [embryoID ': Considering ' label ' curve'] ;
    disp(msg)
    
    % Optimize the placement of the curve7
    lekk = cutCurveIntoLeadingTrailingSegments(curv); 
    lexk = lekk{1} ;
    lesk = lekk{2} ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Chisq
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute chisquared(t) for this curve
    [chisq, chisqn, ~, ssr] = chisquareCurves(lexk, LEX', LEY', (LES.^2)', ...
        true, optimize_trans) ;
    % if allow_rotation
    %     guess = [0, 0, 0]; 
    %     opt = fminsearch(@(x0)curveDifferenceRotTrans(x0, cvfit, xyd), guess, options) ;
    % else
    %     guess = [0, 0]; 
    %     opt = fminsearch(@(x0)curveDifferenceTrans(x0, cvfit, xyd), guess, options) ;
    % end
    
    % Fit chisq to parabola
    % Convert into match and uncertainty, using at least 4 and at most 10
    % timepoints to fit the parabolic minimum.
    chisqn = fillmissing(chisqn, 'movmedian', 5) ;
    if hands_on
        [matchtime, matchtime_unc, fit_coefs, nidx, pidx] = ...
            chisqMinUncInteractiveDomain(chisqn, 4, 10, ssr) ;
    else         
        [tmatch, unc, fit_coefs] = chisqMinUncertainty(chisqn, 4, 10) ;
    end
    matchtime_minutes = matchtime * dt ;
    matchtime_unc_minutes = matchtime_unc * dt ;

    % Plot the result from ChiSq unc
    timedense = 1:0.1:length(chisqn) ;
    % minimimum of the fit is cstar
    cstar = fit_coefs.ystar ;
    close all ;
    fig = figure('visible', 'off') ;
    plot(chisqn, '.-')
    hold on;
    yy = polyval(fit_coefs.p, timedense, fit_coefs.S, fit_coefs.mu) ;
    plot(timedense, yy, '--')
    errorbar(matchtime, cstar, matchtime_unc, 'horizontal')
    ylim([0, 30])
    xlims = xlim() ;
    xlim([0 xlims(2)])
    ylabel('\chi^2 / N')
    xlabel('time [min]')
    title([embryoID, ...
        ': a = ' num2str(matchtime_minutes), ...
        ' \pm ', num2str(matchtime_unc_minutes), '  min'])
    
    % Save into all embryoID matching dirs
    for qq=1:length(eDirs4ID.folders)
        edir = eDirs4ID.folders{qq} ;
        disp(['Saving stripe7_chisq_fit.png into ' edir])
        figfn = fullfile(edir, 'stripe7_chisq_fit.png') ;
        saveas(fig, figfn)
    end
    clf
    
    % Save Chisq timing as mat and txt for all matching embryoID dirs
    for qq=1:length(eDirs4ID.folders)
        edir = eDirs4ID.folders{qq} ;
        disp(['Saving timematch_curve7_chisq.mat into ' edir])
        
        fnmat = fullfile(edir, 'timematch_curve7_chisq.mat') ;
        fntxt = fullfile(edir, 'timematch_curve7_chisq.txt') ;
        save(fnmat, 'matchtime', 'matchtime_unc', 'matchtime_minutes', 'matchtime_unc_minutes')    

        disp(['Saving matchtime to ', fntxt])
        dlmwrite(fntxt, [matchtime_minutes, matchtime_unc_minutes])
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SSR 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Filter to get minimum of SSR
    % windowWidth = 3;  % for quick smoothing
    % kernel = ones(windowWidth, 1) / windowWidth;
    % ssdsm = filter(kernel, 1, ssds);
    [~, minID] = min(ssr) ;
    optID = max(1, minID-6):min(minID+6, length(ssr)) ;
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
    matchtime_minutes = matchtime * dt ;
    matchtime_unc_minutes = matchtime_unc * dt ;
    
    % save an image of the five-pt fit
    close all
    fig = figure('visible', 'off') ;
    % errorbar(trange(optID) * dt / 60, pentad, pentad_unc)
    hold on;
    plot(optID * dt, pentad, 's')
    ttt = min(optID * dt):0.1:max(optID * dt) ;
    tttm = ttt * dt ; % convert to minutes
    plot(tttm, polyval(p, ttt, S), '-')
    plot(matchtime_minutes, polyval(p, matchtime, S), 'o', 'color', green)
    ploterr(matchtime_minutes, polyval(p, matchtime, S), matchtime_unc_minutes, [])
    ylabel('Sum of squared distances')
    xlabel('time [min]')
    title([embryoID, ...
        ': a = ' num2str(matchtime_minutes), ...
        ' \pm ', num2str(matchtime_unc_minutes), '  min'])
    ylims = ylim;
    plot((1:length(ssr))*dt, ssr, 'k.--')
    ylim(ylims)
    fn = fullfile(embryoDir, 'stripe7_ssr_fit.png') ;
    disp(['Saving figure to ' fn])
    saveas(fig, fn)
    close all
    
    % Save timing as mat and txt
    for qq=1:length(eDirs4ID.folders)
        edir = eDirs4ID.folders{qq} ;
        disp(['Saving timematch_curve7_chisq.mat into ' edir])
        fnmat = fullfile(edir, 'timematch_curve7_ssr.mat') ;
        fntxt = fullfile(edir, 'timematch_curve7_ssr.txt') ;
        save(fnmat, 'matchtime', 'matchtime_unc', 'matchtime_minutes', 'matchtime_unc_minutes')

        disp(['Saving matchtime to ', fntxt])
        dlmwrite(fntxt, [matchtime_minutes, matchtime_unc_minutes])
    end
    
end