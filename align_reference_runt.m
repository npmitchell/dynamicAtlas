% align_reference_runt.m 
% Script for aligning dynamic runt nanobody data against each other
% 
% NPMitchell 2020

%% ADD PATHS
gitDir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/' ;
gutDir = fullfile(gitDir, 'gut_matlab') ;
basicsDir = fullfile(gutDir, 'basics') ;
tiffDir = fullfile(gutDir, 'tiff_handling') ;
plottingDir = fullfile(gutDir, 'plotting') ;
codeDir = '/Users/npmitchell/Box/Flies/code/time_alignment_2020/time_align_embryos/' ;
fmDir = fullfile(gutDir, 'toolbox_fast_marching/toolbox_fast_marching/') ;
fmDir2 = fullfile(fmDir, 'mex') ;
fmDir3 = fullfile(fmDir, 'toolbox') ;
addpath(basicsDir) ;
addpath(tiffDir) ;
addpath(plottingDir) ;
addpath(codeDir) ;
addpath(fmDir) ;
addpath(fmDir2) ;
addpath(fmDir3) ;

%% OPTIONS
% Save each runt nanobody (curated) MIP as ./date/cylinder1_max.tif
runtNBodyDir = './Runt-Nanobody/' ;
mipfn = 'cylinder1_max.tif' ;
outdir = './alignment' ;
ssfactor = 4 ;              % subsampling factor before computing corr
% Correlation options
corr_method = 'realspace' ; % realspace or phase method for correlation 
                            % If realspace, does not tranlate, phase allows
                            % dx,dy then computes realspace corr on shifted
                            % image.
stripe7corr_method = 'dist' ; 
corrOutDir = fullfile(outdir, [corr_method '_corr']) ;

dirs2make = {outdir, corrOutDir} ;
for ii = 1:length(dirs2make)
    dir2make = dirs2make{ii} ;
    if ~exist(dir2make, 'dir')
        mkdir(dir2make)
    end
end

% General options
preview = false ;           % display intermediate results
overwrite = false ;         % overwrite previous results
thres = 0.5 ;

% define naming
mipfnBase = mipfn(1:end-4) ;
extn = [sprintf('_ss%02d', ssfactor) '_' corr_method] ;

exptnames = {'202001141730', '202001141943', '202001142033', ...
    '202001150004', '202001210000', '202001210044'} ;
expts = cell(length(exptnames), 1) ;
for ii = 1:length(exptnames)
    expts{ii} = fullfile(runtNBodyDir, exptnames{ii}) ;
end
% expts = subdirs(runtNBodyDir) ;
substr = sprintf('_ss%02d', ssfactor) ;
ssmipfn = [mipfnBase substr '.tif'];

%% Aesthetics
colors = define_colors() ;
blue = colors(1, :) ;
red = colors(2, :) ;
yellow = colors(3, :) ;
purple = colors(4, :) ;
green = colors(5, :) ;
sky = colors(6, :) ;
colors = [yellow; sky] ;

%% Subsample the data
for ii = 1:length(expts)
    disp(['dataset ii = ', num2str(ii)])
    % Load time sequence MIPs of dataset ii
    outfn = fullfile(expts{ii}, ssmipfn) ;
    if ~exist(outfn, 'file') || overwrite
        disp('--> DOWNSAMPLING')
        disp(['Loading MIPs tiff: ' fullfile(expts{ii}, mipfn)])
        mipI = loadtiff(fullfile(expts{ii}, mipfn)) ;
        for qq = 1:size(mipI, 3) 
            % Save each subsampled slice
            imI = imadjustn(squeeze(mipI(:, :, qq))) ;   
            imI = imresize( imI, 1./double(ssfactor), 'bilinear' );
            % write the first page as overwrite mode, then append
            if qq == 1
                imwrite(imI, outfn, 'tiff', 'Compression','none');
            else
                imwrite(imI, outfn, 'tiff', 'Compression','none',...
                    'WriteMode','append');    
            end
        end
    end
end

%% Train on stripe7 for each subsampled data in ilastik

%% Extract leading edge of stripe from Loaded probabilities
minsz = 5e3 ;
maxsz = 1e6 ;
for ii = 1:length(expts)
    disp(['Extracting stripe 7 for expt ' num2str(ii)])
    try
        assert(length(fns) == 1)
    catch
        error('More than  one Probabilities.h5 file found!')
    end
    % Try to load curve
    curvfn = fullfile(expts{ii}, 'stripe7curves.mat') ;
    
    % Check if all curves are present
    if exist(curvfn, 'file') 
        continue_curves = false ;
        load(curvfn, 'stripe7curves')
        for tt = 1:length(stripe7curves)
            if isempty(stripe7curves{tt})
                continue_curves = true ;
            end
        end
    else
        continue_curves = true ;
    end
        
    if continue_curves || overwrite
        % Announce what is going on
        if exist(curvfn, 'file')
            disp('Overwriting curv on disk')
        else
            disp('Curv not found on disk, defining')
        end
        
        if save_images
            % Load the full resolution tiff images for image overlays
            tiffstack_fullres = loadtiff(fullfile(expts{ii}, mipfn)) ;
        end
        
        % Get all timepoints for this dataset
        fns = dir(fullfile(expts{ii}, [ mipfnBase substr '_Probabilities.h5'])) ;
        pfn = fullfile(fns(1).folder, fns(1).name) ;
        disp(['Reading ' pfn])
        dat = h5read(pfn, '/exported_data') ;
        stripeDir = fullfile(expts{ii}, ['stripe7' substr]) ;
        if ~exist(stripeDir, 'dir')
            mkdir(stripeDir)
        end

        % No curve mat file found, so compute the curves from each time image
        curves = cell(size(dat, 4), 1) ;
        
        for tt = 1:size(dat, 4)
            % time string
            timestr = sprintf('%04d', tt) ;
            if mod(tt, 1) == 0
                msg = ['Obtaining cropped data for ' ...
                    timestr '/' num2str(size(dat,4))] ;
                disp(msg)
            end
            
            maskfn = fullfile(stripeDir, ['mask_' timestr '.tif']);  
            if exist(maskfn, 'file')
                bw2 = imread(maskfn);
            else
                disp(['mask does not exist: ' maskfn])
                aux_auto_stripe7
                
                % Check if this automatic way is good enough
                move_on = false ;
                good_enough = false ;
                while ~move_on        
                    % Restart the figure options
                    close all
                    imshow(bw2)
                    msg = [timestr ': Enter=OK, ', ...
                        'Backspace=manual ROI, c=clear, a=anteriorCC, e=erode, d=dilate, r=reuse prev ROI, u=undo'] ;
                    title(msg)
                    disp(msg)
                    button = waitforbuttonpress() ;
                    if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                        move_on = true ;
                        good_enough = true ;
                    elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                        move_on = true;
                        good_enough = false ;
                    elseif button && strcmp(get(gcf, 'CurrentKey'), 'c')
                        % Clear and reset to autogenerate
                        aux_auto_stripe7
                    elseif button && strcmp(get(gcf, 'CurrentKey'), 'r')
                        if exist('BWmask', 'var')
                            % reuse previous mask
                            bw2 = BWmask .* bw ;
                        else
                            disp('No BWmask in memory!')
                        end    
                        imshow(bw2)
                        colorbar()
                        title(msg)
                        disp(msg)
                        button = waitforbuttonpress() ;
                        % good_enough = false ;
                        % if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                        %     disp('Good enough using masking of other centroid objects, saving mask')
                        %     disp(['Saving mask to ' maskfn])
                        %     imwrite(bw2, maskfn) ;
                        %     move_on = true ;
                        % elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                        %     disp('Define manual ROI poly')
                        %     move_on = true ;
                        % end
                    elseif button && (strcmp(get(gcf, 'CurrentKey'), 'e') ||...
                            strcmp(get(gcf, 'CurrentKey'), 'd'))
                        % Enter loop where we continually morphologically
                        % refine by eroding + dilating
                        aux_morphological_op
                        disp('end of routine e/d')
                    elseif button && (strcmp(get(gcf, 'CurrentKey'), 'a') || ...
                            strcmp(get(gcf, 'CurrentKey'), 'p'))
                        aux_anterior_movement
                        disp('end of routine a')
                    elseif button && strcmp(get(gcf, 'CurrentKey'), 'u')
                        % recall current state [ie UNDO]
                        bw = bwprev2 ;
                        disp('end of routine e/d')
                    elseif button
                        disp('Bad button press. Use a,e,d,u,return,backspace')
                    end
                    bwprev2 = bwprev ;
                    bwprev = bw ;
                    disp('outside selection loop')
                end
                close all

                % If not good enough, add manual polygon mask
                if good_enough
                    imwrite(bw2, maskfn)
                else
                    if ~exist(maskfn, 'file')
                        imshow(bw)
                        title(timestr)
                        disp('Select ROI polygon (BWmask)')
                        BWmask = roipoly ;
                        % Check out the result
                        bw2 = BWmask .* bw ;
                        imshow(bw2)
                        title('Close to continue')    
                        waitfor(gcf)
                        
                        disp(['Saving mask to ' maskfn])
                        imwrite(bw2, maskfn) ;                
                    else
                        bw2 = imread(maskfn);
                    end
                end
            end

            % Create the curve by finding the mean true pixel in bw image
            % curvc = zeros(size(bw2, 2), 1) ;
            % for jj = 1:size(bw2, 2)
            %     row = bw2(:, jj) ;
            %     try
            %         curvc(jj) = mean(find(row)) ; % , 1, 'last') ;
            %     catch
            %         disp('no true values in this row. Skipping...')
            %     end
            % end
            % inds = find(curvc) ;
            % curv1 = [curvc(inds) + addx, inds + addy] ;
            % curv2 = [flipud(curvc(inds)) + addx, inds + addy] ;
            % curv = [curv1; curv2] ;
            % [yvals, ind] = sort(curv(:, 2)) ;
            % curv = curv(ind, :) ;
            % plot(curv(:, 1), curv(:, 2), '.')
            % axis equal
            % pause(0.001)

            % Create the curve by finding the exterior boundaries
            % Filter out small regions
            bw2 = bwareafilt(logical(bw2), [minsz maxsz]) ;
            curv = bwboundaries(bw2, 'noholes');
            assert(length(curv) == 1)
            curv = curv{1} ;
            curv = [curv(:, 1) + addx, curv(:, 2) + addy] ;
            
            % Store this curve
            disp(['Adding curve for t=', num2str(tt)])
            curves{tt} = curv ;

            % Save fancy image
            figfn = fullfile(stripeDir, ['stripe7' substr '_' timestr '.png']) ;
            if save_images && ~exist(figfn, 'file')
                im = tiffstack_fullres(:, :, tt) ;
                % Adjust intensity
                % [f,x] = ecdf(im(:));
                % f1 = find(f>cdf_min, 1, 'first');
                % f2 = find(f<cdf_max, 1, 'last');
                % lim = [x(f1) x(f2)];
                % im = mat2gray(double(im),double(lim));
                im = imadjustn(im) ;

                % Now make the figure
                fig = figure('visible', 'off') ;
                imshow(im)
                hold on;
                plot(curv(:, 1)*ssfactor, curv(:, 2)*ssfactor,...
                    '.', 'color', yellow)
                title('Identification of stripe 7')
                saveas(fig, figfn)
            end
            
            if mod(tt, 10) == 0
                % Save the curves as a mat file
                stripe7curves = curves ;
                save(curvfn, 'stripe7curves')
            end
        end

        % Save the curves as a mat file
        stripe7curves = curves ;
        save(curvfn, 'stripe7curves')
    end
end


%% Align one timeseries against another using xcorr2fft
for ii = 1:length(expts)
    disp(['dataset ii = ', num2str(ii)])
    % Load time sequence MIPs of dataset ii
    disp(['Loading MIPs tiff: ' fullfile(expts{ii}, mipfn)])
    mipI = loadtiff(fullfile(expts{ii}, mipfn)) ;
    ntpI = size(mipI, 3) ;
    for jj = 1:length(expts)
        % Define the correlation matrix filename
        ijstr = [ '_%02d_%02d' extn ] ;
        cfn = fullfile(corrOutDir, sprintf(['corr' ijstr '.mat'], ii, jj)) ;
        disp(['Seeking cfn = ' cfn])
        
        % Decide to compute the correlations or not
        if ii < jj && (~exist(cfn, 'file') || overwrite)
            if ~exist(cfn, 'file')
                disp('Computing correlations')
            else
                disp('Overwriting correlations')
            end
            disp(['  dataset jj = ', num2str(jj)])
            % Load time sequence MIPs of dataset jj
            disp(['  Loading MIPs tiff: ' fullfile(expts{jj}, mipfn)])
            mipJ = loadtiff(fullfile(expts{jj}, mipfn)) ;
            ntpJ = size(mipJ, 3) ;

            % Compare dataset ii with all timepoints of dataset ii and jj
            % Store phase correlations in cij, phase correlation offsets
            cij = zeros(ntpI, ntpJ) ;
            dxij = zeros(ntpI, ntpJ) ;
            dyij = zeros(ntpI, ntpJ) ;
            for ti = 1:ntpI
                disp(['  > ti = ' num2str(ti)])
                % grab the timepoint ti of MIPs dataset ii
                imI = imadjustn(squeeze(mipI(:, :, ti))) ;   
                imI = imresize( imI, 1./double(ssfactor), 'bilinear' );
                for tj = 1:ntpJ
                    % if mod(tj, 50) == 0
                    %     disp(['  >> tj = ' num2str(tj)])
                    % end
                    % grab the timepoint tj of MIPs dataset jj
                    imJ = imadjustn(squeeze(mipJ(:, :, tj))) ;
                    imJ = imresize( imJ, 1./double(ssfactor), 'bilinear' );
                    
                    if strcmp(corr_method, 'phase')
                        [dxij(ti, tj), dyij(ti,tj), ~] = ...
                            xcorr2fft(imI, imJ) ;
                        [imIshift, imJshift] = shiftImagesX(dxij(ti, tj), imI, imJ, 0) ;
                        cij(ti, tj) = corr2(imIshift, imJshift) ;
                    else
                        cij(ti, tj) = corr2(imI, imJ) ;
                    end
                    
                    % Preview progress in correlations
                    if mod(tj, ntpJ) == 0 && mod(jj, 5) == 0 && preview 
                        subplot(2,2, 1)
                        imagesc(imI)
                        subplot(2,2, 2)
                        % grab middle imJ
                        imJ = imadjustn(squeeze(mipJ(:, :, round(ntpJ*0.5)))) ;
                        imJ = imresize( imJ, 1./double(ssfactor), 'bilinear' );
                        imagesc(imJ)
                        subplot(2,2,[3,4])
                        plot(1:ntpJ, cij(ti, :)); 
                        hold on
                        title(['corr = ' num2str(cij(ti, tj))])
                        pause(0.1)
                    end
                end
            end
            
            % Save cij, dxij, dyij        
            if strcmp(corr_method, 'phase')
                save(cfn, 'dxij', 'dyij', 'cij') 
            else
                save(cfn, 'cij') 
            end
            
            % Plot heatmap of phase correlations, dxs, and dys
            close all
            imagesc(cij)
            xlabel(['timepoint for dataset ' num2str(jj)])
            ylabel(['timepoint for dataset ' num2str(ii)])
            axis equal
            cfn = fullfile(outdir, sprintf(['cij_' corr_method ijstr '.png'], ii, jj)) ;
            cb = colorbar() ;
            ylabel(cb, 'correlation')
            set(gca,'YDir','normal')
            axis tight
            if strcmp(corr_method, 'phase')
                title('phase correlation intensity')
            else
                title('cross correlation')
            end
            saveas(gcf, cfn)

            if strcmp(corr_method, 'phase')
                % delta x
                clf
                imagesc(dxij)
                xlabel(['timepoint for dataset ' num2str(jj)])
                ylabel(['timepoint for dataset ' num2str(ii)])
                axis equal
                cfn = fullfile(outdir, sprintf(['dxij_' corr_method ijstr '.png'], ii, jj)) ;
                cb = colorbar() ;
                ylabel(cb, '\deltax')
                set(gca,'YDir','normal')
                axis tight
                title('phase correlation offset \deltax')
                saveas(gcf, cfn)        
                % delta y 
                clf
                imagesc(dyij)
                xlabel(['timepoint for dataset ' num2str(jj)])
                ylabel(['timepoint for dataset ' num2str(ii)])
                axis equal
                cfn = fullfile(outdir, sprintf(['dyij_' corr_method ijstr '.png'], ii, jj)) ;
                cb = colorbar() ;
                ylabel(cb, '\deltay')
                set(gca,'YDir','normal')
                axis tight
                title('phase correlation offset \deltay')
                saveas(gcf, cfn)           
            end
        elseif ii > jj 
            disp('Skipping since ii > jj')
        else
            disp('Skipping since already done')
        end
    end
end

%% ID peaks in correlation to draw curves -- possibly use fast marcher?
clearvars imI imJ imA imB
for ii = 1:length(expts)
    disp(['dataset ii = ', num2str(ii)])
    for jj = 1:length(expts)
        % Define the correlation matrix filename
        ijstr = [ '_%02d_%02d' extn ] ;
        cfn = fullfile(corrOutDir, sprintf(['corr' ijstr '.mat'], ii, jj)) ;
        disp(['Seeking cfn = ' cfn])
        cpathfn = fullfile(corrOutDir, sprintf(['cpath' ijstr '.mat'], ii, jj)) ;
        cpathImFn = fullfile(corrOutDir, sprintf(['cpathImFn' ijstr '.mat'], ii, jj)) ;
        
        % Decide to compute the correlations or not
        if ii < jj && (~exist(cpathfn, 'file') || overwrite)
            this_ijstr = sprintf(ijstr, ii, jj) ;

            % Create directory for overlays
            overlayDir = fullfile(corrOutDir, ['overlay' this_ijstr]) ;
            if ~exist(overlayDir, 'dir')
                mkdir(overlayDir) 
            end
            outfnBase = fullfile(overlayDir, ['overlay' this_ijstr '_%04d']) ;
                        
            % load the correlations
            load(cfn, 'cij') ;
            % detect which dimension is smaller
            if size(cij, 2) < size(cij, 1) 
                swap = true ;
                cij = cij' ;
                imA = mipJ ;
                imB = mipI ;
            else
                swap = false ;
                imA = mipI ;
                imB = mipJ ;
            end
            
            % find the path
            tpath = zeros(size(cij, 1), 2) ;
            tpath(:, 1) = 1:size(cij, 1) ;
            for tq = 1:size(cij, 1)
                [~, ind] = max(cij(tq, :)) ;
                tpath(tq, 2) = ind ;
            end
            
            % Use correlation heatmap
            close all
            imagesc(cij);
            hold on;
            plot(tpath(:, 2), tpath(:, 1), 'o') 
            axis equal
            axis tight
            msg = 'Does path look ok? Enter=yes, backspace=no' ;
            xlabel(['time, dataset ', num2str(ii)])
            ylabel(['time, dataset ', num2str(jj)])
            title(msg)
            disp(msg)
            good_button = false ;
            while ~good_button
                button = waitforbuttonpress() ;
                if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                    move_on = true ;
                    good_button = true ;
                elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                    move_on = false;
                    good_button = true ;
                end
            end
            
            % Shortest path 
            Woffset = 0. ;
            while ~move_on
                close all
                imagesc(cij);
                hold on;
                plot(tpath(:, 2), tpath(:, 1), 'ro') 
                plot(tpath(1, 2), tpath(1, 1), 'ks')
                plot(tpath(end, 2), tpath(end, 1), 'k^')
                msg = 'Do endpoints look ok? Enter=yes, backspace=no/select' ;
                title(msg)
                disp(msg)
                xlabel(['time, dataset ', num2str(ii)])
                ylabel(['time, dataset ', num2str(jj)])
                button = waitforbuttonpress() ;
                if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                    move_on = true ;
                    startpt = tpath(1, :) ;
                    endpt = tpath(end, :) ;
                elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                    move_on = false;
                    msg = 'Select startpoint' ;
                    % This function is from Gabriel Peyre for picking pts
                    msg = 'Select endpoint' ;
                    title(msg)
                    disp(msg)
                    [startpt, endpt] = pick_start_end_point(cij) ;
                end
                
                % get shortest path via fastmarching by Gabriel Peyre
                options.propagation_type = 'normal';
                options.Tmax = sum(size(cij))*1.2;
                options.start_points = startpt';
                clf;
                disp('Performing FM');
                options.reduc_factor = 1.0 ;
                options.weight = 1.0;
                %   'W' is the weight matrix (the highest, the slowest the front will move).
                %   'start_points' is a 2 x k array, start_points(:,i) is the ith starting point .
                %   'end_points' is a 2 x 1 array, it is the goal.
                % Define speed of pixels' movement for marching
                % Blur the image a bit to straighten lines, then find curve
                % WW = imgaussfilt((max(cij(:)) - cij) - min(cij(:)), 1) ;
                % WW = imgaussfilt(cij) ;
                % !!!
                WW = imgaussfilt(cij - min(cij(:)) + Woffset) ;
                % WW = WW - min(WW(:)) + 1e-2 ;
                % WW = cij - min(cij(:)) + 1e-2;
                % Compute distance transform
                % didn't get slow version to work
                % [DD,S] = perform_front_propagation_2d_slow(WW',...
                %     [startpt(2);startpt(1)], [endpt(2); endpt(1)], 4000, []);
                %
                [DD,S] = perform_fast_marching(WW', [startpt(2)+1; startpt(1)+1], options) ;
                imagesc(DD)
                button = waitforbuttonpress() ;
                
                % perform_fmstar_2d(WW', spt, ept, options);    
                disp('Extracting Paths');
                stepsize = 0.1 ;
                str_options = [stepsize 10000];
                
                % path extraction
                options.str_options = str_options ;
                options.trim_path = true ;
                options.startpt = [startpt(2), startpt(1)] ;
                options.thres_dist = 2 ;
                str_options = options.str_options ;

                % grad = compute_grad(DD);
                % grad = -perform_vf_normalization(grad);
                % Dx = squeeze(grad(:, :, 1)) ;
                % Dy = squeeze(grad(:, :, 2)) ;

                % figure;
                % subplot(1, 2, 1)
                % imagesc(squeeze(grad(:, :, 1))); title('grad(1)')
                % colormap(bwr)
                % caxis([-1,1])
                % subplot(1, 2, 2)
                % imagesc(Dy); title('Dy')
                % colormap(bwr)
                % caxis([-0.01,0.01])
                % 
                % figure;
                % subplot(1, 2, 1)
                % imagesc(squeeze(grad(:, :, 2))); title('grad(2)')
                % colormap(bwr)
                % caxis([-1,1])
                % subplot(1, 2, 2)
                % imagesc(Dx); title('Dx')
                % colormap(bwr)
                % caxis([-1,1])


                [Dy, Dx] = gradient(DD) ;
                Dx = - Dx ;
                Dy = - Dy ;
                
                normalization = (Dx.^2 + Dy.^2);
                stationary = find(normalization < 1e-6) ;
                normalization(stationary) = 1;
                Dx = Dx .* (1./sqrt(normalization)) ;
                Dy = Dy .* (1./sqrt(normalization)) ;
                
                imagesc(DD'); hold on;
                quiver(Dx', Dy', 1)
                title('DT with gradient')
                button = waitforbuttonpress() ;
                
                % Plot with streamlines
                xx2=xx(1:10:size(xx, 1), 1:10:size(xx, 2));
                yy2=yy(1:10:size(yy, 1), 1:10:size(yy, 2));
                cxyz = stream2(xx, yy, Dx', Dy', xx2', yy2') ;
                imagesc(cij); hold on;
                quiver(xx, yy, Dy', Dx', 0)
                title('cij with gradient')
                plot(startpt(2), startpt(1), 'ko')
                plot(endpt(2), endpt(1), 'ko')
                streamline(cxyz)
                figure(5);
                imagesc(cij); hold on;
                quiver(Dx', Dy', 0)
                plot(startpt(2), startpt(1), 'ko')
                plot(endpt(2), endpt(1), 'ko')
                streamline(cxyz)
                button = waitforbuttonpress() ;

                % path extraction
                % grad = cat(3, Dx, Dy) ;
                % grad = perform_vf_normalization(grad) ;
                % Dx = squeeze(grad(:, :, 1)) ;
                % Dy = squeeze(grad(:, :, 2)) ;
                [xx, yy ] = meshgrid(1:size(Dx, 1), 1:size(Dx, 2)) ;
                % works but is uphill
                % path = stream2(Dx, Dy, endpt(1),endpt(2), str_options);
                % messing around here
                path = stream2(xx, yy, Dx', Dy', endpt(2)-1, endpt(1)-1, str_options);
                % works but does not connect to startpt
                % path = stream2(Dy, Dx, endpt(1), endpt(2), str_options) ;
                % path = stream2(-Dx', -Dy', startpt(1), startpt(2), str_options) ;
                path = path{1} ;
                path = [path(:, 2), path(:, 1)] ;
               
                figure(1)
                imagesc(DD') ;
                title('DT')
                colorbar()
                hold on;
                plot(path(:, 2), path(:, 1), '.-')
                plot(startpt(:, 2), startpt(:, 1), 'o')
                plot(endpt(:, 2), endpt(:, 1), 'o')
                button = waitforbuttonpress() ;

                figure(2);
                imagesc(Dx') ;
                colormap(bwr)
                title('Dx')
                colorbar()
                %caxis([-1,1])
                hold on;
                plot(path(:, 2), path(:, 1), '.-')
                plot(startpt(:, 2), startpt(:, 1), 'o')
                plot(endpt(:, 2), endpt(:, 1), 'o')
                button = waitforbuttonpress() ;

                figure(3);
                imagesc(Dy') ;
                title('Dy')
                colormap(bwr)
                colorbar()
                %caxis([-1,1])
                hold on; 
                plot(path(:, 2), path(:, 1), '.-')
                plot(startpt(:, 2), startpt(:, 1), 'o')
                plot(endpt(:, 2), endpt(:, 1), 'o')
                button = waitforbuttonpress() ;

                figure(4);
                imagesc(cij) ;
                title('cij')
                colormap(bwr)
                colorbar()
                caxis([-1,1])
                hold on; 
                plot(path(:, 2), path(:, 1), '.-')
                plot(startpt(:, 2), startpt(:, 1), 'o')
                plot(endpt(:, 2), endpt(:, 1), 'o')
                button = waitforbuttonpress() ;

                % TEST
                xx2=xx(1:10:size(xx, 1), 1:10:size(xx, 2));
                yy2=yy(1:10:size(yy, 1), 1:10:size(yy, 2));
                cxyz = stream2(xx, yy, Dx', Dy', xx2', yy2') ;
                figure(5);
                imagesc(cij); hold on;
                quiver(Dx', Dy', 0)
                streamline(cxyz)
                button = waitforbuttonpress() ;
                                
                % cpath = compute_geodesic(DD, flipud(ept), opt);
                % Clip the path near the start point
                d2start = vecnorm(path-startpt,2,2) ;
                path = path(d2start > thres_dist, :) ;
                [tmp, ind] = unique(path(:, 1)) ;
                inds = sort(ind) ;
                path = path(inds, :) ;
                cpath = [endpt; path; startpt] ;
               
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Visualize resulting path
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                close all
                imagesc(cij);
                hold on;
                plot(startpt(2), startpt(1), 'ro') 
                plot(endpt(2), endpt(1), 'ks')
                plot(cpath(:, 2), cpath(:, 1), '-') ;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Convert cpath to tpath
                % Remember that the first dimension (rows) is shorter
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                ntimeptsA = round(endpt(1) - startpt(1) + 1) ; 
                xxx = round(startpt(1)):round(endpt(1)) ;
                yyy = interp1(flipud(cpath(:, 1)), ...
                    flipud(cpath(:, 2)), xxx) ;
                tpath = round([xxx; yyy]') ;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Continue visualization
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                plot(tpath(:, 2), tpath(:, 1), 'o') ;
                xlabel(['time, dataset ', num2str(ii)])
                ylabel(['time, dataset ', num2str(jj)])
                title('Path ok? Enter=yes, backspace=no/redo')
                button = waitforbuttonpress() ;
                if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                    move_on = true ;
                elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                    move_on = false;
                    disp('starting over with path detection')
                    Woffset = input('Set weight offset = ') ;
                    disp(' --> new Woffset')
                    if isempty(Woffset)
                        Woffset = input('Set weight offset = ') ;
                    end
                end
            end
            
            % Save image
            figfn = fullfile(corrOutDir, ['correspondence' this_ijstr '.png']) ;
            disp(['Saving ' figfn])
            title(['Correspondences between ' num2str(ii) ' and ' num2str(jj) ])
            saveas(gcf, figfn) ;
            
            % For each match, plot overlay as cyan magenta
            overlayDir = fullfile(corrOutDir, ['overlay' this_ijstr]) ;
            if ~exist(overlayDir, 'dir')
                mkdir(overlayDir) 
            end
            outfn = fullfile(overlayDir, ['overlay' this_ijstr '_%04d.png']) ;
            
            % Save overlays
            close all
            for qq = 1:size(cij, 1)
                if mod(qq, 10) == 0
                    disp(['Saving overlaid stripes ' num2str(qq)])
                end
                set(gcf, 'visible', 'off')
                ima = double(imadjustn(squeeze(imA(:, :, qq)))) ;
                imb = double(imadjustn(squeeze(imB(:, :, tpath(qq, 2))))) ;
                ima = uint8(255*mat2gray(ima, [0 max(ima(:))])) ;
                imb = uint8(255*mat2gray(imb, [0 max(imb(:))])) ;
                [shiftdat, shiftfixed] = shiftImagesX(0, ima, imb, 0) ;
                rgb = zeros(size(ima, 1), size(ima, 2), 3, 'uint8') ;
                rgb(:, :, 1) = ima ;
                rgb(:, :, 2) = imb ;
                rgb(:, :, 3) = min(ima + imb, 255) ;
                % save rgb image
                imshow(rgb)
                title(['Correspondence between ' num2str(ii) ' and ' num2str(jj)])
                set(gcf, 'PaperUnits', 'centimeters');
                set(gcf, 'PaperSize', [16, 16]) ;
                saveas(gcf, sprintf(outfn, qq))
            end
            
            % Save the path
            time_correspondences = tpath;
            save(cpathfn, 'time_correspondences') ;
        end
    end
end


%% Align Stripe 7 for each frame -- correlations in position
for ii = 1:length(expts)
    disp(['dataset ii = ', num2str(ii)])
    % Load time sequence MIPs of dataset ii
    curvIfn = fullfile(expts{ii}, 'stripe7curves.mat') ;
    disp(['  Loading curves: ' curvIfn])
    curvI = load(curvIfn, 'stripe7curves');
    curvI = curvI.stripe7curves ;
    ntpI = length(curvI) ;
    
    for jj = 1:length(expts)
        % Define the correlation matrix filename
        ijstr = [ '_%02d_%02d' extn ] ;
        cfn = fullfile(outdir, sprintf(['corr_stripe7' ijstr '.mat'], ii, jj)) ;
        disp(['Seeking cfn = ' cfn])
        
        % Decide to compute the correlations or not
        if ii < jj && (~exist(cfn, 'file') || overwrite)
            if ~exist(cfn, 'file')
                disp('Computing stripe7 correlations')
            else
                disp('Overwriting stripe7 correlations')
            end
            disp(['  dataset jj = ', num2str(jj)])
            % Load time sequence MIPs of dataset jj
            curvJfn = fullfile(expts{jj}, 'stripe7curves.mat') ;
            disp(['  Loading curves: ' curvJfn])
            curvJ = load(curvJfn, 'stripe7curves');
            curvJ = curvJ.stripe7curves ;
            ntpJ = length(curvJ) ;

            % Compare dataset ii with all timepoints of dataset ii and jj
            % Store phase correlations in cij, phase correlation offsets
            cij = zeros(ntpI, ntpJ) ;
            dxij = zeros(ntpI, ntpJ) ;
            dyij = zeros(ntpI, ntpJ) ;
            for ti = 1:ntpI
                disp(['  > ti = ' num2str(ti)])
                % grab the timepoint ti of MIPs dataset ii
                cI = curvI{ti} ;   
                for tj = 1:ntpJ
                    % if mod(tj, 50) == 0
                    %     disp(['  >> tj = ' num2str(tj)])
                    % end
                    % grab the timepoint tj of MIPs dataset jj
                    cJ = curvJ{tj} ;
                    if strcmp(stripe7corr_method, 'dist')
                        d1 = pdist2(cI, cJ);
                        mindist1 = min(d1, [], 2);
                        d2 = pdist2(cJ, cI);
                        mindist2 = min(d2, [], 2);
                        cij(ti, tj) = mean(sqrt(mindist1)) * mean(sqrt(mindist2)) ;
                    else
                        cij = correlation ;
                    end
                    
                    % Preview progress in correlations
                    if mod(tj, ntpJ) == 0 && mod(jj, 5) == 0 && preview 
                        subplot(2,2, 1)
                        imagesc(imI)
                        subplot(2,2, 2)
                        % grab middle imJ
                        imJ = imadjustn(squeeze(mipJ(:, :, round(ntpJ*0.5)))) ;
                        imJ = imresize( imJ, 1./double(ssfactor), 'bilinear' );
                        imagesc(imJ)
                        subplot(2,2,[3,4])
                        plot(1:ntpJ, cij(ti, :)); 
                        hold on
                        title(['corr = ' num2str(cij(ti, tj))])
                        pause(0.1)
                    end
                end
            end
            
            % Save cij, dxij, dyij        
            save(cfn, 'cij') 
            
            % Plot heatmap of phase correlations, dxs, and dys
            close all
            imagesc(cij)
            xlabel(['timepoint for dataset ' num2str(jj)])
            ylabel(['timepoint for dataset ' num2str(ii)])
            axis equal
            cfn = fullfile(outdir, ...
                sprintf(['stripe7cij_' stripe7corr_method ijstr '.png'], ii, jj)) ;
            cb = colorbar() ;
            cmap = colormap ;
            colormap(flipud(cmap))
            if strcmp(stripe7corr_method, 'dist')
                ylabel(cb, '$\langle d_{ij} \rangle$', 'Interpreter', 'Latex')
                title('Runt stripe 7 $\langle d_{ij} \rangle$', 'Interpreter', 'Latex')
            else
                ylabel(cb, 'correlation')
                title('Runt stripe 7 correlation')
            end
            set(gca,'YDir','normal')
            axis tight
            saveas(gcf, cfn)
        elseif ii > jj 
            disp('Skipping since ii > jj')
        else
            disp('Skipping since already done')
        end
        
    end
end
disp('done')


%% Minimize difference between timing points
% {tj} = {dtj} + frame, where 
% Fit line of best fit to each ridge extraction

