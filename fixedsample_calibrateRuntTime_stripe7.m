%% Determine time of runt image from depth of last stripe
% example fn: /mnt/crunch/Atlas_Data/final_for_box/WT/Runt/201908021230
%
% Use Runt to time each image against Runt Nanobody reference data
save_fancy = true; 
allow_rotation = false ;
overwrite = false ;
overwrite_ROI = false ;
prepend = 'Max_Cyl1_2_000001_c1_rot_scaled_view1.tif' ;
cdf_min = 0.01 ;
cdf_max = 0.999 ;
sigma = 5 ;
kernel = 5*sigma;
step = 1 ;
exten = sprintf('_sigma%03d_step%03d', sigma, step) ;
thres = 0.77 ;
width = 500 ;
dt = 75 ;
manualROI = false ;
wexten = ['_width' num2str(width)] ;

% Add paths
addpath('./polyparci')
addpath('./ploterr')
gutDir = '/data/code/gut_matlab/' ;
addpath(gutDir)
addpath('/data/code/time_align_embryos/lookup/')

%% Build the lookuptable
a = lookupTable ;
a = a. ;
lut = a.map('Runt') ;
lut2 = a.map('Slp') ;

%% Plotting
[colors, names] = define_colors ;
yellow = colors(3, :) ;
green = colors(5, :) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make smoothed images if they don't exist
% First convert Runt to smoothed image at const luminosity
label = 'Runt' ;
for kk = 1:length(lut.folders)
    embryoDir = lut.folders{kk} ;
    sfn = fullfile(embryoDir, [label '_smooth' exten '.tif']) ;
    
    if ~exist(sfn, 'file')
        % Build the dir of this experiment
        todo = dir(fullfile(embryoDir, [prepend '2' postpend])) ;
        if length(todo) > 1
            error('here')
        end
        im = imread(fullfile(todo(ii).folder, todo(ii).name)) ;
        
        % Adjust intensity
        [f,x] = ecdf(im(:));
        f1 = find(f>cdf_min, 1, 'first');
        f2 = find(f<cdf_max, 1, 'last');
        lim = [x(f1) x(f2)];
        scale = mat2gray(double(im),double(lim));
        smoothim = imfilter(scale, fspecial('gaussian', kernel, sigma));

        % Save imsc
        imwrite(imsc, sfn)
    else
        disp('smoothed image already on disk')
    end
   
end
disp('done')

%% Process probabilities in iLastik -- train on stripe7 versus not a stripe

%% Load dynamic data Eve curves near where Runt should lie
load('../time_calibration_stripe7/eve/stripe7curves_tmin25_tmax48.mat', 'stripe7curves') ;
evecurvs = stripe7curves ;
trange = 25:41 ;

%% Load probabilities
for kk = 1:length(lut.folders)
    
    embryoDir = lut.folders{kk} ;
    embryoID = lut.embryoIDs{kk} ;
    
    stripefn = fullfile(embryoDir, ['Runt_stripe7curve' wexten '.mat']) ;
    if ~exist(stripefn, 'file') || overwrite

        pfn = fullfile(embryoDir, [label '_smooth' exten '_Probabilities.h5']) ;
        dat = h5read(pfn, '/exported_data') ;
        midx = round(0.5 * size(dat, 2)) ;
        midy = round(0.5 * size(dat, 3)) ;
        dcrop = squeeze(dat(2, midx:end, midy-width:midy+width)) ;
        addx = midx ;
        addy = midy - width ;

        % Extract the leading curve
        maskfn = fullfile(embryoDir, ['Runt_mask' wexten '.tif']);  
        if exist(maskfn, 'file') && ~overwrite_ROI
            % Apply manual ROI to thresholded image with default thres
            BWmask = imread(maskfn) ;

            figure ;
            % Start from scratch with default thres
            thres0 = thres ;
            move_on = false ;
            while ~move_on
                datbw = false(size(dcrop)) ;
                datbw(dcrop > thres0) = true ;
                bw = bwareafilt(datbw, 8) ;
                bw = bwareafilt(bw, [1e3 1e8]) ;
                imshow(bw .* BWmask)
                title(['Loaded BWmask. Enter to continue, <-> for threshold change. thres = ', num2str(thres0)])
                button = waitforbuttonpress() ;
                if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                    move_on = true ;
                elseif button && strcmp(get(gcf, 'CurrentKey'), 'rightarrow')
                    thres0 = min(1, thres0 + 0.05) ;
                elseif button && strcmp(get(gcf, 'CurrentKey'), 'leftarrow')
                    thres0 = max(0, thres0 - 0.05) ;
                end
            end
            
            % Check out the result
            bw = BWmask .* bw ;
            % Detect and adjust the curve
            curv = zeros(size(bw, 2), 1) ;
            for jj = 1:size(bw, 2)
                row = bw(:, jj) ;
                try
                    curv(jj) = find(row, 1, 'last') ;
                catch
                    disp('no true values in this row. Skipping...')
                end
            end
            inds = find(curv) ;
            if any(isnan(curv))
                inds = inds(~isnan(curv)) ;
            end
            curv_adj = [curv(inds) + addx, inds + addy] ;

            % Prepare the plot
            close all
            imshow(1-dcrop)
            hold on
            plot(curv, 'o-')
            curv
            title(['Identification of stripe 7: ' embryoID '. <Enter> to save'])
            % Check if the result is good
            button = waitforbuttonpress() ;
            if button && ~strcmp(get(gcf, 'CurrentKey'), 'return')
                error('Exiting since non-Return keypress in figure check')
            else
                disp('Saving figure and stripe...')
            end
        else
            thres0 = thres ;
            move_on = false ;
            while ~move_on
                datbw = false(size(dcrop)) ;
                datbw(dcrop > thres0) = true ;
                posterior = bwareafilt(datbw, 1) ;
                datbw(posterior) = false ;
                bw = bwareafilt(datbw, 8) ;
                bw = bwareafilt(bw, [1e4 1e8]) ;
                cc = bwconncomp(bw) ;
                rp = regionprops(cc) ;
                centry = zeros(length(rp), 1) ;
                for qq = 1:length(rp)
                    centry(qq) = rp(qq).Centroid(2) ;
                end
                [~, ind] = max(centry) ;
                obw = false(size(dcrop)) ;
                if isempty(cc.PixelIdxList)
                    if thres0 ~= thres
                        disp("No pixels detected, returning to original threshold")
                        thres0 = thres ;
                    else
                        disp('No pixels with default thres. Changing thres to random value')
                        thres0 = rand(1) ;
                    end
                else
                    obw(cc.PixelIdxList{ind}) = true ;

                    % Detect and adjust the curve
                    curv = zeros(size(obw, 2), 1) ;
                    for jj = 1:size(obw, 2)
                        row = obw(:, jj) ;
                        try
                            curv(jj) = find(row, 1, 'last') ;
                        catch
                        end
                    end
                    inds = find(curv) ;
                    if any(isnan(curv))
                        inds = inds(~isnan(curv)) ;
                    end
                    curv_adj = [curv(inds) + addx, inds + addy] ;

                    % Save an image of the identification
                    close all
                    fig = figure ;
                    imshow(1-dcrop)
                    hold on
                    plot(curv, 'o-')
                    title(['Enter to continue. <-> to adjust threshold. Backspace for ROI. thres=' num2str(thres0)])
                    button = waitforbuttonpress() ;
                    if button && strcmp(get(fig, 'CurrentKey'), 'return')
                        move_on = true ;
                    elseif button && strcmp(get(fig, 'CurrentKey'), 'rightarrow')
                        thres0 = min(1, thres0 + 0.05) ;
                    elseif button && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
                        thres0 = max(0, thres0 - 0.05) ;
                    elseif button && strcmp(get(fig, 'CurrentKey'), 'backspace')
                        % Do manual ROI
                        ofn = fullfile(embryoDir, [label '_smooth' exten '.tif']) ;
                        oim = imread(ofn) ;
                        figure ;
                        imshow(oim) 
                        
                        figure ;
                        % Start from scratch with default thres
                        thres0 = thres ;
                        % There used to be an option to mask the posterior
                        % midgut, but this is difficult to use in practice
                        % mask_posterior = true ;
                        while ~move_on
                            datbw = false(size(dcrop)) ;
                            datbw(dcrop > thres0) = true ;
                            % posterior = bwareafilt(datbw, 1) ;
                            % if mask_posterior
                            %     datbw(posterior) = false ;
                            % end
                            bw = bwareafilt(datbw, 8) ;
                            bw = bwareafilt(bw, [1e3 1e8]) ;
                            imshow(bw)
                            title(['Enter to continue, <-> for threshold change. thres = ', num2str(thres0)])
                            button = waitforbuttonpress() ;
                            if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                                move_on = true ;
                            elseif button && strcmp(get(gcf, 'CurrentKey'), 'rightarrow')
                                thres0 = min(1, thres0 + 0.05) ;
                            elseif button && strcmp(get(gcf, 'CurrentKey'), 'leftarrow')
                                thres0 = max(0, thres0 - 0.05) ;
                            end
                        end
                        % Select the ROI
                        close all
                        imshow(bw)
                        title(['Select ROI: ' embryoID '\n' embryoDir])
                        BWmask = roipoly ;
                        disp(['Saving mask to ' maskfn])
                        imwrite(BWmask, maskfn) ;

                        % Check out the result
                        bw = BWmask .* bw ;
                        imshow(bw)
                        title('Close to continue')    

                        % Detect and adjust the curve
                        curv = zeros(size(bw, 2), 1) ;
                        for jj = 1:size(bw, 2)
                            row = bw(:, jj) ;
                            try
                                curv(jj) = find(row, 1, 'last') ;
                            catch
                                disp('no true values in this row. Skipping...')
                            end
                        end
                        inds = find(curv > 0) ;
                        if any(isnan(curv))
                            inds = inds(~isnan(curv)) ;
                        end
                        curv_adj = [curv(inds) + addx, inds + addy] ;

                        % Prepare the plot
                        close all
                        imshow(1-dcrop)
                        hold on
                        plot(curv, 'o-')
                    end
                end
            end
            % Save the image
            title(['Identification of stripe 7: ' embryoID '. <Enter> to save'])
            % Check if the result is good
            if button && ~strcmp(get(gcf, 'CurrentKey'), 'return')
                error('Exiting since non-Return keypress in figure check')
            else
                disp('Saving figure and stripe')
            end
            saveas(gcf, fullfile('../time_calibration_stripe7/', [embryoID '.png']))
            close all
        end

        % depth = max(curv) - min(curv) ;
        stripe7curve = curv_adj; 
        save(stripefn, 'stripe7curve')
    
        % Save fancy image
        if save_fancy
            % Load image
            im1 = imread(fullfile(ch1(1).folder, ch1(1).name)) ;

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

