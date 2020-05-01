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

%% Build the lookuptable
a = lookupTable ;
a = a.buildLookup('') ;
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
chanl = 'Runt' ;
for kk = 1:length(lut.folders)
    exptDir = lut.folders{kk} ;
    sfn = fullfile(exptDir, [channel '_smooth' exten '.tif']) ;
    
    if ~exist(sfn, 'file')
        % Build the dir of this experiment
        todo = dir(fullfile(exptDir, [prepend '2' postpend])) ;
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
    
    exptDir = lut.folders{kk} ;
    exptName = strsplit(lut.folders{kk}, '/') ;
    exptName = exptName{end} ;
    
    stripefn = fullfile(exptDir, ['Runt_stripe7curve' wexten '.mat']) ;
    if ~exist(stripefn, 'file') || overwrite

        pfn = fullfile(exptDir, [channel '_smooth' exten '_Probabilities.h5']) ;
        dat = h5read(pfn, '/exported_data') ;
        midx = round(0.5 * size(dat, 2)) ;
        midy = round(0.5 * size(dat, 3)) ;
        dcrop = squeeze(dat(2, midx:end, midy-width:midy+width)) ;
        addx = midx ;
        addy = midy - width ;

        % Extract the leading curve
        maskfn = fullfile(exptDir, ['Runt_mask' wexten '.tif']);  
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
            title(['Identification of stripe 7: ' exptName '. <Enter> to save'])
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
                        ofn = fullfile(exptDir, [channel '_smooth' exten '.tif']) ;
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
                        title(['Select ROI: ' exptName '\n' exptDir])
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
            title(['Identification of stripe 7: ' exptName '. <Enter> to save'])
            % Check if the result is good
            if button && ~strcmp(get(gcf, 'CurrentKey'), 'return')
                error('Exiting since non-Return keypress in figure check')
            else
                disp('Saving figure and stripe')
            end
            saveas(gcf, fullfile('../time_calibration_stripe7/', [exptName '.png']))
            close all
        end

        % depth = max(curv) - min(curv) ;
        stripe7curve = curv_adj; 
        save(stripefn, 'stripe7curve')

        if strcmp(stripefn, '/Users/npmitchell/Box/Flies/WT/Ftz_Paired_Runt/201905091640/Runt_stripe7curve.mat') || kk == 8
            error('here')
        end
    
        % Save fancy image
        if save_fancy
            cavailable = lut.channelnums{kk} ;

            % Load each channel
            ch1 = dir(fullfile(exptDir, [prepend num2str(cavailable(1)) postpend])) ;
            im1 = imread(fullfile(ch1(1).folder, ch1(1).name)) ;

            % Adjust intensity
            [f,x] = ecdf(im1(:));
            f1 = find(f>cdf_min, 1, 'first');
            f2 = find(f<cdf_max, 1, 'last');
            lim = [x(f1) x(f2)];
            im1 = mat2gray(double(im1),double(lim));

            ch2 = dir(fullfile(exptDir, [prepend num2str(cavailable(2)) postpend])) ;
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
                ch3 = dir(fullfile(exptDir, [prepend num2str(cavailable(2)) postpend])) ;
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
            saveas(fig, fullfile('../time_calibration_stripe7/', ['rgb_' exptName '.png']))

        end
        
    end
end


%% Now match to time domain Eve data
for kk = 1:length(lut.folders)
    exptDir = lut.folders{kk} ;
    exptName = strsplit(lut.folders{kk}, '/') ;
    exptName = exptName{end} ;

    % Which adjusted curve
    load(fullfile(exptDir, ['Runt_stripe7curve' wexten '.mat']), 'stripe7curve') ;
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
            if allow_rotation
                guess = [0, 0, 0]; 
                opt = fminsearch(@(x0)curveDifferenceRotTrans(x0, cvfit, xyd), guess, options) ;
            else
                guess = [0, 0]; 
                opt = fminsearch(@(x0)curveDifferenceTrans(x0, cvfit, xyd), guess, options) ;
            end
            
            % rotate and translate based on optimal params
            if allow_rotation
                th = opt(1) ;
                dxy = [opt(2) opt(3)] ; 
                cvfit = ([ cos(th), -sin(th); sin(th), cos(th) ] * cvfit')' ;
            else
                dxy = [opt(1) opt(2)] ;
            end
            cvfit = cvfit + dxy ;
            % Find sum of squared diffs via interpolation
            [xd, index] = unique(xyd(:, 1)) ;
            ydint = interp1(xd, xyd(index, 2), cvfit(:, 1), 'pchip') ;
            ssds = [ssds sum((cvfit(:, 2) - ydint).^2)] ;
            
            % store timestamp
            tstamps  = [tstamps qq] ;
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
    % plot(tstamps, ssds)
    % plot(ssds)
    % plot(tstamps(optID), ssds(optID))    
    % Fit to y = a x^2 + b x + c
    [p, S] = polyfit(trange(optID) * dt - trange(minID), pentad, 2) ;
    a = p(1) ;
    b = p(2) ;
    % [y_fit,delta] = polyval(p, trange(pentx) * dt, S);
    ci = polyparci(p, S, 0.6827) ;
    a_unc = a - ci(1, 1) ;
    b_unc = b - ci(1, 2) ;
    % Find minimum
    if a > 0
        matchtime = -b / (2 * a) + t0 ;
        dmdb = -1 / (2 * a) ;
        dmda = 0.5 * b / a^2 ; 
        matchtime_unc = sqrt((dmda * a_unc)^2 + (dmdb * b_unc)^2) ;
    else
        % the quadratic fit is so bad that it is upside down
        matchtime = trange(minID) ;
        matchtime_unc = NaN ;
    end
    matchtime_minutes = matchtime / 60 ;
    matchtime_unc_minutes = matchtime_unc / 60 ;
    
    % save an image of the five-pt fit
    close all
    fig = figure('visible', 'off') ;
    % errorbar(trange(optID) * dt / 60, pentad, pentad_unc)
    hold on;
    plot(trange(optID) * dt / 60, pentad, 's')
    ttt = min(trange(optID) * dt):0.1:max(trange(optID) * dt) ;
    tttm = ttt / 60 ; % convert to minutes
    plot(tttm, polyval(p, ttt- t0, S), '-')
    plot(matchtime_minutes, polyval(p, matchtime - t0, S), 'o', 'color', green)
    ploterr(matchtime_minutes, polyval(p, matchtime - t0, S), matchtime_unc_minutes, [])
    ylabel('Sum of squared distances')
    xlabel('time [min]')
    title(['a = ' num2str(matchtime_minutes) ' \pm ' num2str(matchtime_unc_minutes) '  min'])
    ylims = ylim;
    plot(trange * dt / 60, ssds, 'k.--')
    ylim(ylims)
    fn = fullfile(exptDir, [exptName '_stripefit_unc_polyfit.png']) ;
    disp(['Saving figure to ' fn])
    saveas(fig, fn)
    close all
    
    % Save timing as mat and txt
    fnmat = ['timematch_curve7' texten '.mat'] ;
    fntxt = ['timematch_curve7' texten '.txt'] ;
    save(fnmat, 'matchtime', 'matchtime_unc', 'matchtime_minutes', 'matchtime_unc_minutes')
    
    disp(['Saving matchtime to ', fntxt])
    dlmwrite(fntxt, [matchtime_minutes, matchtime_unc_minutes])
    
end

