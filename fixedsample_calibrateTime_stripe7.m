%% Determine time of runt image from depth of last stripe
% example fn: /mnt/crunch/Atlas_Data/final_for_box/WT/Runt/201908021230
%
% Use Runt to time each image against Runt Nanobody reference data
clearvars
clc
save_fancy = true; 
allow_rotation = false ;
overwrite = true ;
overwrite_ROI = false ;

% Options: Which data to analyze
label = 'Runt' ;  % Whether to use Runt-Nanobody or Eve data
genoDir = './WT' ; 
prepend = 'MAX_Cyl1_2_000000_c1_rot_scaled_view1' ;
exten = '.tif' ;
if strcmp(label, 'Runt')
    timerfn = 'timematch_RuntNanobody_stripe7.mat' ;
elseif strcmp(label, 'Eve')
    timerfn = 'timematch_EveYFP_stripe7.mat' ;
else
    error('Label is not Runt or Eve. What stripe7 is this?')
end
cdf_min = 0.01 ;
cdf_max = 0.999 ;
sigma = 5 ;
kernel = 5*sigma;
step = 1 ;
sigmastep = sprintf('_sigma%03d_step%03d', sigma, step) ;
dt = 60 ;

%% Add paths
tlaDir = '/Users/npmitchell/Box/Flies/code/time_alignment_2020/time_align_embryos';
if ~exist(tlaDir, 'dir')
    tlaDir = '/data/code/time_align_embryos/' ;
end
gutDir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/gut_matlab/' ;
if ~exist(gutDir, 'dir')
    gutDir = '/data/code/gut_matlab/mfl2/gut_matlab/' ;
end
addpath(fullfile(tlaDir, 'polyparci'))
addpath(fullfile(tlaDir, 'ploterr'))
addpath(fullfile(tlaDir, 'lookup'))
addpath(fullfile(tlaDir, 'fixed_stripe7'))
addpath(gutDir)
addpath(fullfile(gutDir, 'plotting')) ;

%% Build the lookuptable
a = lookupTable() ;
a = a.buildLookup(genoDir, timerfn, prepend, exten) ;
lut = a.map(label) ;

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
probfn = ['sigma020_step001/smooth/', label, '_smooth', sigmastep,...
    '_bin_Probabilities.h5'] ;
for kk = 1:length(lut.folders)
    
    filename = lut.names{kk} ;
    embryoDir = lut.folders{kk} ;
    embryoID = lut.embryoIDs{kk} ;
    
    stripefn = fullfile(embryoDir, 'Runt_stripe7curve.mat') ;
    if ~exist(stripefn, 'file') || overwrite

        pfn = fullfile(embryoDir, probfn) ;
        disp(['Opening ' pfn])
        dat = h5read(pfn, '/exported_data') ;
        midx = round(0.5 * size(dat, 2)) ;
        midy = round(0.5 * size(dat, 3)) ;
        dcrop = squeeze(dat(2, midx:end, :)) ;
        addx = midx ;
        addy = 0 ;

        % Extract the leading curve
        maskfn = fullfile(embryoDir, ['Runt_mask.tif']);  
        if ~exist(maskfn, 'file') || ~exist(stripefn, 'file') || overwrite
            thres = 0.5 ;
            minmaxsz = [5e3, 1e8] ;
            curv = extractStripeEdges(dat, maskfn, thres, minmaxsz, embryoID) ;
            % depth = max(curv) - min(curv) ;
            stripe7curve = curv ; 
            stripe7curve_frac = curv ;
            stripe7curve_frac(:, 1) = double(stripe7curve_frac(:, 1)) / double(size(dat, 2)) ;
            stripe7curve_frac(:, 2) = double(stripe7curve_frac(:, 2)) / double(size(dat, 3)) ;
            save(stripefn, 'stripe7curve', 'stripe7curve_frac')
        else
            curv = load(stripefn, 'stripe7curve') ;
        end

    
        % Save fancy image
        if save_fancy
            % Load image and all other images of this embryo
            estruct = a.findEmbryo(embryoID) ;
            
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
                        
            % Now make the figure
            fig = figure('visible', 'off') ;
            imshow(combined)
            hold on;
            imwrite(combined, fullfile(embryoDir, [alllabels 'rgb_' embryoID '.png']))
            plot(curv(:, 1), curv(:, 2), 'o-', 'color', stripecolor)
            % plot(curv' + midx, (1:length(curv)) + midy - width, 'o-', 'color', yellow)
            title([embryoID ': Identification of stripe 7'])
            saveas(fig, fullfile(embryoDir, [label 'stripe7_rgb_' embryoID '.png']))
        end        
    end
end

%% Now match to time domain Runt Nanobody or EveYFP data
for kk = 1:length(lut.folders)
    embryoDir = lut.folders{kk} ;
    embryoID = lut.embryoIDs{kk} ;

    % Which adjusted curve
    load(fullfile(embryoDir, stripefn), 'stripe7curve') ;
    curv = stripe7curve ;
    
    % Match to dynamic data curv
    ssds = [] ; 
    tstamps = [] ;
    options = optimset('MaxIter', 25, 'TolX', 1e-2) ; 
    % Other options: ('FunctionTolerance', 1e-7, 'Display','iter','PlotFcns',@optimplotfval);
    % Consider each timepoint
    for qq = trange
        % compare to this curve
        if ~isempty(evecurvs{qq})
            msg = [exptName ': Considering eve curve for timept ' num2str(qq)] ;
            disp(msg)
            
            xyd = evecurvs{qq} ;
            xyd = xyd(~isnan(xyd(:, 1)), :) ;
            xyd = [xyd(:, 2), xyd(:, 1) ];
            cvfit = [curv(:, 2), curv(:, 1)];

            % Optimize the placement of the curve7
            
            
            % Compute chisquared(t) for this curve
            [chisq, chisqn, ssr] = ...
                chisquareCurves(curv, refcurvsX, refcurvsY, refvars, smooth_var) ;
            % if allow_rotation
            %     guess = [0, 0, 0]; 
            %     opt = fminsearch(@(x0)curveDifferenceRotTrans(x0, cvfit, xyd), guess, options) ;
            % else
            %     guess = [0, 0]; 
            %     opt = fminsearch(@(x0)curveDifferenceTrans(x0, cvfit, xyd), guess, options) ;
            % end
            
        end
    end
    
    % Filter to get minimum
    % windowWidth = 3;  % for quick smoothing
    % kernel = ones(windowWidth, 1) / windowWidth;
    % ssdsm = filter(kernel, 1, ssds);
    [~, minID] = min(ssds) ;
    optID = max(1, minID-2):min(minID+2, length(ssds)) ;
    % pentad = ssdsm(optID) ;    
    pentad = ssds(optID) ;    
    t0 = trange(minID) ;
    
    % Fit to parabola
    xx = stripe(:, 1) / wAP ;
    yy = stripe(:, 2) / wDV ;
    segs = cutCurveIntoLeadingTrailingSegments([xx, yy]) ;
    lead = segs{1} ;
    curv = lead(:, [2 1]) ;
    [chisq, chisqn, ssr] = chisquareCurves(curv, refcurvX, refcurvsY, refvars, smooth_var) ;
    % error('here')
    [tmatch, unc, fit_coefs] = chisqMinUncertainty(chisqn, 4, 10) ;

    % Plot the result
    timedense = 1:0.1:length(chisqn) ;
    % minimimum of the fit is cstar
    cstar = fit_coefs.ystar ;
    plot(chisqn, '.-')
    hold on;
    yy = polyval(fit_coefs.p, timedense, fit_coefs.S, fit_coefs.mu) ;
    plot(timedense, yy, '--')
    errorbar(tmatch, cstar, unc, 'horizontal')
    ylim([0, 30])
    xlims = xlim() ;
    xlim([0 xlims(2)])
    ylabel('\chi^2 / N')
    xlabel(['time [min], dataset' num2str(ii)])
    title(['dataset ' num2str(ii) ': $t=$' sprintf('%04d', qq)], ...
        'Interpreter', 'Latex')
    figfn = fullfile(embryoDir, sprintf('chisq_snap_%04d.png', qq)) ;
    saveas(fig, figfn)
    clf
        
    matchtime_minutes = matchtime / 60 ;
    matchtime_unc_minutes = matchtime_unc / 60 ;
    
    % Save timing as mat and txt
    fnmat = ['timematch_curve7' texten '.mat'] ;
    fntxt = ['timematch_curve7' texten '.txt'] ;
    save(fnmat, 'matchtime', 'matchtime_unc', 'matchtime_minutes', 'matchtime_unc_minutes')
    
    disp(['Saving matchtime to ', fntxt])
    dlmwrite(fntxt, [matchtime_minutes, matchtime_unc_minutes])
    
end

