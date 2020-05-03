%% Determine time of runt image from depth of last stripe
% example fn: /mnt/crunch/Atlas_Data/final_for_box/WT/Runt/201908021230
%
% Use Runt to time each image against Runt Nanobody reference data
clearvars
clc
save_fancy = true; 
allow_rotation = false ;
overwrite = false ;
overwrite_ROI = false ;

% Options: Which data to analyze
genoDir = './WT' ; 
prepend = 'MAX_Cyl1_2_000000_c1_rot_scaled_view1' ;
exten = '.tif' ;
timerfn = 'timematch_RuntNanobody_stripe7.mat' ;
cdf_min = 0.01 ;
cdf_max = 0.999 ;
sigma = 5 ;
kernel = 5*sigma;
step = 1 ;
sigmastep = sprintf('_sigma%03d_step%03d', sigma, step) ;
thres = 0.77 ;
width = 500 ;
dt = 75 ;
manualROI = false ;
wexten = ['_width' num2str(width)] ;

%% Add paths
tlaDir = '/Users/npmitchell/Box/Flies/code/time_alignment_2020/time_align_embryos';
addpath(fullfile(tlaDir, 'polyparci'))
addpath(fullfile(tlaDir, 'ploterr'))
gutDir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/gut_matlab/' ;
addpath(gutDir)
addpath(fullfile(gutDir, 'plotting')) ;
addpath(fullfile(tlaDir, 'lookup'))
addpath(fullfile(tlaDir, 'fixed_stripe7'))

%% Build the lookuptable
a = lookupTable ;
a = a.buildLookup(genoDir, timerfn, prepend, exten) ;
lut = a.map('Runt') ;

%% Plotting
[colors, names] = define_colors ;
yellow = colors(3, :) ;
green = colors(5, :) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assumes we have already run EnsembleStripe to make smoothed images 

%% Create overlay of all stripes for each embryo
lut2 = a.map('Even_Skipped') ;
for kk = 1:length(lut.folders)
    lut2.loadSecondar()
end


%% Process probabilities in iLastik -- train on stripe7 versus not a stripe

%% Identify stripes if not done so already
label = 'Runt' ;
for kk = 1:length(lut.folders)
    
    filename = lut.names{kk} ;
    embryoDir = lut.folders{kk} ;
    embryoID = lut.embryoIDs{kk} ;
    
    stripefn = fullfile(embryoDir, ['Runt_stripe7curve.mat']) ;
    if ~exist(stripefn, 'file') || overwrite

        pfn = fullfile(embryoDir, sigmastep(2:end), ...
            [label '_smooth' sigmastep '_Probabilities.h5']) ;
        disp(['Opening ' pfn])
        dat = h5read(pfn, '/exported_data') ;
        midx = round(0.5 * size(dat, 2)) ;
        midy = round(0.5 * size(dat, 3)) ;
        dcrop = squeeze(dat(2, midx:end, :)) ;
        addx = midx ;
        addy = midy - width ;

        % Extract the leading curve
        maskfn = fullfile(embryoDir, ['Runt_mask.tif']);  
        if exist(maskfn, 'file') && ~overwrite_ROI
            curv = extractStripeEdges(dat, maskfn) ;
        end

        % depth = max(curv) - min(curv) ;
        stripe7curve = curv ; 
        save(stripefn, 'stripe7curve')
    
        % Save fancy image
        if save_fancy
            % Load image
            im1 = imread(fullfile(embryoDir, filename)) ;
            
            % Adjust intensity
            [f,x] = ecdf(im1(:));
            f1 = find(f>cdf_min, 1, 'first');
            f2 = find(f<cdf_max, 1, 'last');
            lim = [x(f1) x(f2)];
            im1 = mat2gray(double(im1),double(lim));

            ch2 = dir(fullfile(embryoDir, [prepend num2str(cavailable(2)) postpend])) ;
            im2 = imread(fullfile(ch2(1).folder, ch2(1).name)) ;

            % Adjust intensity
            [f,x] = ecdf(im2(:));
            f1 = find(f>cdf_min, 1, 'first');
            f2 = find(f<cdf_max, 1, 'last');
            lim = [x(f1) x(f2)];
            im2 = mat2gray(double(im2),double(lim));

            combined = zeros(size(im2, 1), size(im2, 2), 3) ;
            combined(:, :, 1) = im1 ;
            combined(:, :, 2) = im2 ;

            if length(cavailable) > 2
                ch3 = dir(fullfile(embryoDir, [prepend num2str(cavailable(2)) postpend])) ;
                im3 = imread(fullfile(ch3(1).folder, ch3(1).name)) ;
                % Adjust intensity
                [f,x] = ecdf(im3(:));
                f1 = find(f>cdf_min, 1, 'first');
                f2 = find(f<cdf_max, 1, 'last');
                lim = [x(f1) x(f2)];
                im3 = mat2gray(double(im3),double(lim));
                combined(:, :, 3) = im3 ;
            end

            % Now make the figure
            fig = figure('visible', 'off') ;
            imshow(combined)
            hold on;
            plot(curv_adj(:, 1), curv_adj(:, 2), 'o-', 'color', yellow)
            % plot(curv' + midx, (1:length(curv)) + midy - width, 'o-', 'color', yellow)
            title('Identification of stripe 7')
            saveas(fig, fullfile('../time_calibration_stripe7/', ['rgb_' embryoID '.png']))
        end        
    end
end

%% Now match to time domain Eve data
for kk = 1:length(lut.folders)
    embryoDir = lut.folders{kk} ;
    embryoID = lut.embryoIDs{kk} ;

    % Which adjusted curve
    load(fullfile(embryoDir, ['Runt_stripe7curve' wexten '.mat']), 'stripe7curve') ;
    curv_adj = stripe7curve ;
    
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
            cvfit = [curv_adj(:, 2), curv_adj(:, 1)];

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

