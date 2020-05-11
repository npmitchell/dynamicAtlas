
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
probfn = [ mipfnBase substr '_Probabilities.h5'] ;
minsz = 5e3 ;
maxsz = 1e6 ;

for ii = 1:length(expts)
    disp(['Extracting stripe 7 for expt ' num2str(ii)])
    fns = dir(fullfile(expts{ii}, probfn)) ;
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
        continue_curves = false ;
        error('here')
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
            
            % Here, we cut the image in half, assuming that the time we are
            % interested in finishes before stripe 7 extends halfway across
            % the dorsal side
            % Find the center pixel (rounded to nearest integer)
            midx = round(0.5 * size(dat, 2)) ;
            % Crop the image to right and center
            dcrop = squeeze(dat(1, midx:end, :, tt)) ;
            dsz = [size(dat, 2), size(dat, 3)] ;
            addx = midx ;
            addy = 0 ;

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
                    if exist('bwprev', 'var')
                        bwprev2 = bwprev ;
                    end
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


%% Extract cross correlations between experiments
corrDatOutDir = fullfile(corrOutDir, 'crosscorrelation') ;
if ~exist(corrDatOutDir, 'dir')
    mkdir(corrDatOutDir)
end

% consider each experiment pair
for ii = 1:length(expts)
    disp(['dataset ii = ', num2str(ii)])
    
    loaded_ii = false ;
    for jj = 1:length(expts)
        % Define the correlation matrix filename
        ijstr = [ '_%02d_%02d' extn ] ;
        cfn = fullfile(corrDatOutDir, sprintf(['corr' ijstr '.mat'], ii, jj)) ;
        disp(['Seeking cfn = ' cfn])
        
        % Decide to compute the correlations or not 
        if ~exist(cfn, 'file') || overwrite
            
            % Load dataset ii now since it is needed
            if ~loaded_ii 
                % Load time sequence MIPs of dataset ii
                disp(['Loading MIPs tiff: ' fullfile(expts{ii}, mipfn)])
                mipI = loadtiff(fullfile(expts{ii}, mipfn)) ;
                ntpI = size(mipI, 3) ;
                loaded_ii = true ;
            end
            
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
            cfn = fullfile(corrDatOutDir, sprintf(['cij' ijstr '.png'], ii, jj)) ;
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
                cfn = fullfile(corrDatOutDir, sprintf(['dxij_' corr_method ijstr '.png'], ii, jj)) ;
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
                cfn = fullfile(corrDatOutDir, sprintf(['dyij_' corr_method ijstr '.png'], ii, jj)) ;
                cb = colorbar() ;
                ylabel(cb, '\deltay')
                set(gca,'YDir','normal')
                axis tight
                title('phase correlation offset \deltay')
                saveas(gcf, cfn)           
            end
        else
            disp('Skipping since already done')
        end
    end
end
disp('Done with cross-correlation')

%% ID peaks in correlation to draw curves using fast marching
clear imI imJ imA imB
corrImOutDir = fullfile(corrOutDir, 'correspondence_images') ;
corrPathOutDir = fullfile(corrOutDir, 'cpath') ;
dirs2make = {corrImOutDir, corrPathOutDir} ;
for di = 1:length(dirs2make) 
    dir2make = dirs2make{di} ;
    if ~exist(dir2make, 'dir')
        mkdir(dir2make)
    end
end

for ii = 1:length(expts)
    disp(['dataset ii = ', num2str(ii)])
    
    loaded_ii = false ;
    
    for jj = 1:length(expts)
        % Define the correlation matrix filename
        ijstr = [ '_%02d_%02d' extn ] ;
        cfn = fullfile(corrDatOutDir, sprintf(['corr' ijstr '.mat'], ii, jj)) ;
        disp(['Seeking cfn = ' cfn])
        cpathfn = fullfile(corrPathOutDir, sprintf(['cpath' ijstr '.mat'], ii, jj)) ;
        
                
        % Decide to compute the correlations or not
        if ii == jj && (~exist(cpathfn, 'file') || overwrite)
            
            % Load dataset ii if not already done
            if ~loaded_ii
                % Load time sequence MIPs of dataset ii
                disp(['Loading MIPs tiff: ' fullfile(expts{ii}, mipfn)])
                mipI = loadtiff(fullfile(expts{ii}, mipfn)) ;
                loaded_ii = true ;
            end
            
            % Self-correspondence
            this_ijstr = sprintf(ijstr, ii, jj) ;
            tpath = [(1:size(mipI, 3))', (1:size(mipI, 3))'] ;
            time_correspondences = tpath;
            save(cpathfn, 'time_correspondences') ;
        elseif ~exist(cpathfn, 'file') || overwrite
            
            % Load dataset ii if not already done
            if ~loaded_ii
                % Load time sequence MIPs of dataset ii
                disp(['Loading MIPs tiff: ' fullfile(expts{ii}, mipfn)])
                mipI = loadtiff(fullfile(expts{ii}, mipfn)) ;
                loaded_ii = true ;
            end
            
            % Load time sequence MIPs of dataset jj
            disp([' --> Loading MIPs tiff for jj=', num2str(jj), ...
                ': ', fullfile(expts{jj}, mipfn)])
            mipJ = loadtiff(fullfile(expts{jj}, mipfn)) ;
            this_ijstr = sprintf(ijstr, ii, jj) ;

            % Create directory for overlays
            overlayDir = fullfile(corrOutDir, ['overlay' this_ijstr]) ;
            if ~exist(overlayDir, 'dir')
                mkdir(overlayDir) 
            end
            outfnBase = fullfile(overlayDir, ['overlay' this_ijstr ]) ;
                        
            % load the correlations
            disp(['Loading cfn from: ' cfn])
            load(cfn, 'cij') ;
            % % detect which dimension is smaller
            % if size(cij, 2) < size(cij, 1) 
            %     swap = true ;
            %     cij = cij' ;
            % else
            %     swap = false ;
            % end
            
            % find the path (rough path via maxima) 
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
            ntpguess = min(size(cij)) ;
            xtmp = (1:ntpguess) + tpath(1, 2) ;
            ytmp = (1:ntpguess) + tpath(1, 1) ;
            plot(xtmp, ytmp, 'k--')
            xtmp = tpath(end, 2) - (0:ntpguess-1) ;
            ytmp = tpath(end, 1) - (0:ntpguess-1) ;
            plot(xtmp, ytmp, 'k--')
            axis equal
            axis tight
            msg = 'Does path look ok? Enter=yes, Backspace=no, n/Delete=No correspondence' ;
            xlabel(['time, dataset ', num2str(ii)])
            ylabel(['time, dataset ', num2str(jj)])
            title(msg)
            disp(msg)
            good_button = false ;
            abort = false ;
            while ~good_button
                button = waitforbuttonpress() ;
                if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                    move_on = true ;
                    good_button = true ;
                elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                    move_on = false;
                    good_button = true ;
                    % Save this 
                    title(['Initial correspondence between ' num2str(ii) ' and ' num2str(jj)])
                    saveas(gcf, fullfile(corrImOutDir, ['correspondence' this_ijstr '_initial.png']))
                elseif button && (strcmp(get(gcf, 'CurrentKey'), 'delete') || strcmp(get(gcf, 'CurrentKey'), 'n') )
                    abort = true ;
                    good_button = true ;
                    move_on = true
                end
            end
            
            % Shortest path 
            Woffset = 0. ;
            cij_exponent = 1. ; 
            while ~move_on
                close all
                imagesc(cij);
                hold on;
                plot(tpath(:, 2), tpath(:, 1), 'ro') 
                plot(tpath(1, 2), tpath(1, 1), 'ks')
                plot(tpath(end, 2), tpath(end, 1), 'k^')
                xtmp = (1:ntpguess) + tpath(1, 2) ;
                ytmp = (1:ntpguess) + tpath(1, 1) ;
                plot(xtmp, ytmp, 'k--')
                xtmp = tpath(end, 2) - (0:ntpguess-1) ;
                ytmp = tpath(end, 1) - (0:ntpguess-1) ;
                plot(xtmp, ytmp, 'k--')
                msg = 'Do endpoints look ok? Enter=yes, Backspace=no/select' ;
                title(msg)
                disp(msg)
                xlabel(['time, dataset ', num2str(ii)])
                ylabel(['time, dataset ', num2str(jj)])
                startpt = tpath(1, :) ;
                % Guess start/endpt to be first/last tpath points
                endpt = tpath(end, :) ;
                button = waitforbuttonpress() ;
                if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                    move_on = true ;
                    disp('set startpt/endpt:')
                    startpt
                    endpt
                elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                    move_on = false;
                    msg = 'Start timept is ok? Enter=yes, Backspace=no' ;
                    scatter(startpt(2), startpt(1), 100, 'r', 'filled')
                    disp(msg)
                    title(msg)
                    button = waitforbuttonpress() ;
                    if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                        disp('Great, how about end point?')
                    elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                        % New guess for startpoint is maximum along column
                        [~, startpt(1, 1)] = max(cij(:, 1)) ;
                        scatter(startpt(2), startpt(1), 100, 'r', 'filled')
                        msg = 'Ok, how about startpt now? Enter=yes, Backspace=no' ;
                        disp(msg)
                        title(msg)
                        button = waitforbuttonpress() ;
                        if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                            disp('Great, how about end point?')
                        elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                            msg = 'Find the startpoint manually by clicking' ;
                            disp(msg)
                            title(msg)
                            startpt = ginput(1) ;
                            startpt = [startpt(2) startpt(1)] ;
                        end
                    end
                    msg = 'End timept is ok? Enter=yes, Backspace=no' ;
                    scatter(startpt(2), startpt(1), 100, 'g', 'filled')
                    scatter(endpt(2), endpt(1), 100, 'r', 'filled')
                    disp(msg)
                    title(msg)
                    button = waitforbuttonpress() ;
                    if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                        disp('Great, all done.')
                    elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                        % New guess for startpoint is maximum along column
                        [~, endpt(1, 1)] = max(cij(:, end)) ;
                        endpt(1, 2) = size(cij, 2) ;
                        scatter(endpt(2), endpt(1), 100, 'r', 'filled')
                        msg = 'Ok, how about startpt now? Enter=yes, Backspace=no' ;
                        disp(msg)
                        title(msg)
                        button = waitforbuttonpress() ;
                        if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                            disp('Great, all done.')
                        elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                            endpt = ginput(1) ;
                            endpt = [endpt(2) endpt(1)] ;
                        end
                    end
                else
                    error('bad button press')
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
                WW = imgaussfilt((cij - min(cij(:))).^ cij_exponent + Woffset) ;
                % WW = WW - min(WW(:)) + 1e-2 ;
                % WW = cij - min(cij(:)) + 1e-2;
                % Compute distance transform
                % didn't get slow version to work
                % [DD,S] = perform_front_propagation_2d_slow(WW',...
                %     [startpt(2);startpt(1)], [endpt(2); endpt(1)], 4000, []);
                %
                [DD,S] = perform_fast_marching(WW', [startpt(2)+0.1; startpt(1)+0.1], options) ;
                
                % Check DD
                % imagesc(DD)
                % button = waitforbuttonpress() ;
                
                % perform_fmstar_2d(WW', spt, ept, options);    
                disp('Extracting Paths');
                stepsize = 0.1 ;
                thres_dist = 2 ;
                str_options = [stepsize 10000];
                
                % path extraction
                options.str_options = str_options ;
                options.trim_path = true ;
                options.startpt = [startpt(2), startpt(1)] ;
                options.thres_dist = thres_dist ;
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


                % Compute gradient
                [Dy, Dx] = gradient(DD) ;
                Dx = - Dx ;
                Dy = - Dy ;
                normalization = (Dx.^2 + Dy.^2);
                stationary = find(normalization < 1e-6) ;
                normalization(stationary) = 1;
                Dx = Dx .* (1./sqrt(normalization)) ;
                Dy = Dy .* (1./sqrt(normalization)) ;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                disp('Show gradient')
                imagesc(DD'); hold on;
                quiver(Dx', Dy', 1)
                title('DT with gradient')
                ccaxis = caxis ; % get current caxis
                caxis([0, min(3, ccaxis(2))])
                cb = colorbar() ;
                ylabel(cb, 'DT') ;
                axis equal 
                axis tight
                pause(0.001)
                saveas(gcf, [outfnBase '_DD.png'])
                
                % Plot with streamlines
                disp('Saving cijstream')
                [xx, yy ] = meshgrid(1:size(Dx, 1), 1:size(Dx, 2)) ;
                xx2=xx(1:10:size(xx, 1), 1:10:size(xx, 2));
                yy2=yy(1:10:size(yy, 1), 1:10:size(yy, 2));
                cxyz = stream2(xx, yy, Dx', Dy', xx2', yy2') ;
                imagesc(cij); hold on;
                quiver(xx, yy, Dy', Dx', 0)
                title('cij with gradient')
                plot(startpt(2), startpt(1), 'ko')
                plot(endpt(2), endpt(1), 'ko')
                title(['Cross-correlation for datasets ' num2str(ii) ' and ' num2str(jj)])
                streamline(cxyz)
                axis equal 
                axis tight
                pause(0.001)
                saveas(gcf, [outfnBase '_cijstream.png'])

                % path extraction
                % grad = cat(3, Dx, Dy) ;
                % grad = perform_vf_normalization(grad) ;
                % Dx = squeeze(grad(:, :, 1)) ;
                % Dy = squeeze(grad(:, :, 2)) ;
                % works but is uphill
                % path = stream2(Dx, Dy, endpt(1),endpt(2), str_options);
                % messing around here
                path = stream2(xx, yy, Dx', Dy', endpt(2)-1, endpt(1)-1, str_options);
                % works but does not connect to startpt
                % path = stream2(Dy, Dx, endpt(1), endpt(2), str_options) ;
                % path = stream2(-Dx', -Dy', startpt(1), startpt(2), str_options) ;
                path = path{1} ;
                path = [path(:, 2), path(:, 1)] ;
               
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                disp('Saving DT')
                figure(1)
                imagesc(DD') ;
                title('DT')
                colorbar()
                hold on;
                plot(path(:, 2), path(:, 1), '.-', 'color', yellow)
                plot(startpt(:, 2), startpt(:, 1), 'ko')
                plot(endpt(:, 2), endpt(:, 1), 'ko')
                xlabel(['time, dataset ', num2str(ii)])
                ylabel(['time, dataset ', num2str(jj)])
                title(['DT for datasets ' num2str(ii) ' and ' num2str(jj)])
                axis equal 
                axis tight
                pause(0.001)
                saveas(gcf, [outfnBase '_DT.png'])

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                disp('Saving Dx')
                figure(2);
                imagesc(Dx') ;
                colormap(bwr)
                title('Dx')
                colorbar()
                %caxis([-1,1])
                hold on;
                plot(path(:, 2), path(:, 1), '.-', 'color', yellow)
                plot(startpt(:, 2), startpt(:, 1), 'ko')
                plot(endpt(:, 2), endpt(:, 1), 'ko')
                xlabel(['time, dataset ', num2str(ii)])
                ylabel(['time, dataset ', num2str(jj)])
                title(['\partial_xW for datasets ' num2str(ii) ' and ' num2str(jj)])
                axis equal 
                axis tight
                pause(0.001)
                saveas(gcf, [outfnBase '_Dx.png'])

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                disp('Saving Dy')
                figure(3);
                imagesc(Dy') ;
                title('Dy')
                colormap(bwr)
                colorbar()
                %caxis([-1,1])
                hold on; 
                plot(path(:, 2), path(:, 1), '.-', 'color', yellow)
                plot(startpt(:, 2), startpt(:, 1), 'ko')
                plot(endpt(:, 2), endpt(:, 1), 'ko')
                xlabel(['time, dataset ', num2str(ii)])
                ylabel(['time, dataset ', num2str(jj)])
                title(['\partial_yW for datasets ' num2str(ii) ' and ' num2str(jj)])
                axis equal 
                axis tight
                pause(0.001)
                saveas(gcf, [outfnBase '_Dy.png'])

                disp('Saving cij')
                figure(4);
                imagesc(cij) ;
                title('cij')
                colormap(bwr)
                colorbar()
                caxis([-1,1])
                hold on; 
                plot(path(:, 2), path(:, 1), '.-', 'color', yellow)
                plot(startpt(:, 2), startpt(:, 1), 'ko')
                plot(endpt(:, 2), endpt(:, 1), 'ko')
                xlabel(['time, dataset ', num2str(ii)])
                ylabel(['time, dataset ', num2str(jj)])
                title(['Correlation between dataset ' num2str(ii) ' and ' num2str(jj)])
                axis equal 
                axis tight
                pause(0.001)
                saveas(gcf, [outfnBase '_cij.png'])

                % TEST
                disp('Saving streamlines')
                xx2=xx(1:10:size(xx, 1), 1:10:size(xx, 2));
                yy2=yy(1:10:size(yy, 1), 1:10:size(yy, 2));
                cxyz = stream2(xx, yy, Dx', Dy', xx2', yy2') ;
                plot(startpt(:, 2), startpt(:, 1), 'ko')
                plot(endpt(:, 2), endpt(:, 1), 'ko')
                figure(5);
                imagesc(cij); hold on;
                quiver(Dx', Dy', 0)
                streamline(cxyz)
                xlabel(['time, dataset ', num2str(ii)])
                ylabel(['time, dataset ', num2str(jj)])
                title(['Correlation between dataset ' num2str(ii) ' and ' num2str(jj)])
                pause(1)
                saveas(gcf, [outfnBase '_streamlines.png'])
                                
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
                    flipud(cpath(:, 2)), xxx, 'pchip') ;
                tpath = [xxx; yyy]' ;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Continue visualization
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                plot(tpath(:, 2), tpath(:, 1), 'o') ;
                xlabel(['time, dataset ', num2str(ii)])
                ylabel(['time, dataset ', num2str(jj)])
                cb = colorbar() ;
                ylabel(cb, 'correlation')            
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
                    
                    cij_exponent = input('Set cij potential exponent = ') ;
                    disp(' --> new cij exponent')
                    if isempty(cij_exponent)
                        cij_exponent = input('Set cij potential exponent = ') ;
                    end
                end
            end
            
            if ~abort
                % % Swap back if we swapped I<->J. Do this before plotting
                % if swap
                %     tpath = [tpath(:, 2), tpath(:, 1)];
                % end
                
                % Save Woffset and cij_exponent        
                cpath_param_fn = fullfile(corrPathOutDir, sprintf(['cpath_params' ijstr '.mat'], ii, jj)) ;
                save(cpath_param_fn, 'Woffset', 'cij_exponent')
                
                % Save image
                figfn = fullfile(corrImOutDir, ['correspondence' this_ijstr '.png']) ;
                disp(['Saving ' figfn])
                title(['Correspondences between ' num2str(ii) ' and ' num2str(jj) ])
                axis equal
                axis tight
                saveas(gcf, figfn) ;

                % For each match, plot overlay as cyan magenta
                overlayDir = fullfile(corrPathOutDir, ['overlay' this_ijstr]) ;
                if ~exist(overlayDir, 'dir')
                    mkdir(overlayDir) 
                end
                outfn = fullfile(overlayDir, ['overlay' this_ijstr '_%04d.png']) ;

                % Save overlays
                close all
                disp('Saving overlays...')
                for qq = 1:length(tpath) 
                    if ~exist(sprintf(outfn, qq), 'file')
                        if mod(qq, 10) == 0
                            disp(['Saving overlaid stripes ' num2str(qq)])
                        end
                        set(gcf, 'visible', 'off')  
                        ima = double(imadjustn(squeeze(mipI(:, :, max(1, min(size(tpath, 1), round(tpath(qq, 1)))))))) ;
                        imb = double(imadjustn(squeeze(mipJ(:, :, max(1, min(size(tpath, 1), round(tpath(qq, 2)))))))) ;
                        ima = uint8(255*mat2gray(ima, [0 max(ima(:))])) ;
                        imb = uint8(255*mat2gray(imb, [0 max(imb(:))])) ;
                        [shiftdat, shiftfixed] = shiftImagesX(0, ima, imb, 0) ;
                        rgb = zeros(size(ima, 1), size(ima, 2), 3, 'uint8') ;
                        rgb(:, :, 1) = ima ;
                        rgb(:, :, 2) = imb ;
                        rgb(:, :, 3) = min(ima + imb, 255) ;
                        % save rgb image
                        imshow(rgb)
                        titlestr = ['Correspondence between ' num2str(ii) ' and ' num2str(jj)] ;
                        titlestr = [titlestr ', $W$=' num2str(Woffset) ': $t_i=$' num2str(tpath(qq,1))] ;
                        title(titlestr, 'Interpreter', 'Latex')
                        axis equal
                        axis tight
                        set(gcf, 'PaperUnits', 'centimeters');
                        set(gcf, 'PaperSize', [16, 16]) ;
                        saveas(gcf, sprintf(outfn, qq))
                    end
                end

                % Save the path
                time_correspondences = tpath;
                save(cpathfn, 'time_correspondences') ;
            end
        end
    end
end
disp('Done with correspondence curves')


%% Align Stripe 7 for each frame -- correlations in position
stripe7CorrDatDir = fullfile(outdir, sprintf('stripe7corr_ss%02d', ssfactor)) ;
if ~exist(stripe7CorrDatDir, 'dir')
    mkdir(stripe7CorrDatDir)
end

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
        ijstr = [ '_%02d_%02d_' stripe7corr_method ] ;
        cfn = fullfile(stripe7CorrDatDir, sprintf(['corr_stripe7' ijstr '.mat'], ii, jj)) ;
        disp(['Seeking cfn = ' cfn])
        
        % Decide to compute the correlations or not
        if ~exist(cfn, 'file') || overwrite
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
                if mod(ti, 20) == 0 
                    disp(['  > ti = ' num2str(ti)])
                end
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
            cfn = fullfile(stripe7CorrDatDir, ...
                sprintf(['stripe7cij' ijstr '.png'], ii, jj)) ;
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
disp('done with correspondence curves')


%% Minimize difference between timing points
% {tj} = {dtj} + frame
% Fit line of best fit to each ridge extraction using time-time-corr cell

timelineDir = fullfile(outdir, sprintf('timeline_ss%02d_realspacecorr', ssfactor)) ;
if ~exist(timelineDir, 'dir')
    mkdir(timelineDir)
end

% Build ttc: time-time correlation cell
ttcfn = fullfile(timelineDir, 'ttc.mat') ;
if exist(ttcfn, 'file') && ~overwrite
    load(ttcfn, 'ttc') ;
else
    ttc = cell(length(expts), 1) ;
    for ii = 1:length(expts)
        % create the cell
        ttc{ii} = cell(length(expts), 1) ;
    end
    % Populate the correspondences into cell from mat files
    for ii = 1:length(expts)
        disp(['dataset ii = ', num2str(ii)])
        for jj = ii:length(expts)
            % Define the correlation matrix filename
            ijstr = [ '_%02d_%02d' extn ] ;
            cfn = fullfile(corrOutDir, sprintf(['corr' ijstr '.mat'], ii, jj)) ;
            disp(['Seeking cfn = ' cfn])
            cpathfn = fullfile(corrPathOutDir, sprintf(['cpath' ijstr '.mat'], ii, jj)) ;
            if exist(cpathfn, 'file')
                load(cpathfn, 'time_correspondences');
                ttc{ii}{jj} = time_correspondences ;
            end
        end
    end
    
    % Save the time-time corr cell
    save(ttcfn, 'ttc')
end

% Plot the dynamics
close all
for use_offset = [true false]
    offset = 0 ;
    for ii = 1:length(ttc)
        subplot(round(length(ttc)*0.5)+1, 2, ii)
        % Do earlier dataset correspondences
        for jj = 1:ii-1
            if ~isempty(ttc{jj}{ii})
                % displace vertically to match
                if use_offset
                    i0 = ttc{jj}{ii}(1, 2) ;
                    j0 = ttc{jj}{ii}(1, 1) ;
                    offset = i0 - j0 ;
                end
                plot(ttc{jj}{ii}(:, 2), ...
                    ttc{jj}{ii}(:, 1) + offset,...
                    '.-', 'color', colorset(jj, :))
                hold on;
            end
        end

        % Do this dataset to itself
        for jj = ii 
            plot(ttc{ii}{jj}(:, 1), ttc{ii}{jj}(:, 2), '.', ...
                'color', colorset(jj, :))
            hold on;
        end

        % Do later dataset correspondences
        for jj = ii+1:length(ttc)
            if ~isempty(ttc{ii}{jj})
                % displace vertically to match
                if use_offset
                    i0 = ttc{ii}{jj}(1, 1) ;
                    j0 = ttc{ii}{jj}(1, 2) ;
                    offset = i0 - j0 ;
                end
                plot(ttc{ii}{jj}(:, 1), ...
                    ttc{ii}{jj}(:, 2) + offset,...
                    '.-', 'color', colorset(jj, :))
                hold on;
            end
        end

        % labels
        ylabel(['$\tau(t_' num2str(ii) ')$' ], 'Interpreter', 'Latex')
        if ii == 5 || ii == 6
            xlabel('time, $t_i$', 'Interpreter', 'latex')
        end
        xlim([1, 150])
        ylim([1, 150]) 
    end

    % Legend
    subplot(round(length(ttc)*0.5)+1, 2, [length(ttc) + 1, length(ttc)+2])
    for ii = 1:length(ttc)
        plot([-1],[-1],'.-', 'color', colorset(ii, :))
        hold on
    end
    legend({'1','2','3','4','5','6'}, 'orientation', 'horizontal')
    xlim([0, 1])
    ylim([0, 1])
    axis off
    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [0 0 16 30]) ;
    version_of_matlab = version ;
    if use_offset
        figurefn = fullfile(timelineDir, 'time_correspondences_offset.png') ;
    else
        figurefn = fullfile(timelineDir, 'time_correspondences.png') ;
    end
    if contains(version_of_matlab, '2020')
        exportgraphics(gcf, figurefn)
    else
        saveas(gcf, figurefn)
    end

    close all
end

%% Relax timepoints to reference curve (time of dataset #hard)
% so-called 'hard reference' is the master timeline
to_relax = 1:length(ttc) ;
to_relax = setdiff(to_relax, hard) ;

% Fit other curves to the hard reference timeline
% Tau0 is a map from timeline ti to reference timeline
% to use, do yy = ppval(pp, xq);
for ii = 1:length(ttc)
    if ii <= hard
        pp = spline(ttc{ii}{hard}(:, 1), ttc{ii}{hard}(:, 2)) ;
    elseif ii > hard
        pp = spline(ttc{hard}{ii}(:, 2), ttc{hard}{ii}(:, 1)) ;
    end
    if ii == 1
        tau0 = [pp] ;
    else
        tau0 = [tau0, pp] ;
    end
end

%% Build i_tau0j master timeline lookup table

% Check for existing result
itau0jfn = fullfile(timelineDir, 'i_tau0j.mat') ;
if exist(itau0jfn, 'file') && ~overwrite
    load(itau0jfn, 'i_tau0j', 'ntp_tot', 'addt', 'ntps')
else
    % First get #tps in each dataset, which is stored in diagonals of ttc
    ntps = zeros(length(ttc), 1) ;
    addt = zeros(length(ttc), 1) ;
    tpid = [] ;
    for ii = 1:length(ttc)
        addt(ii) = sum(ntps) ;
        ntps(ii) = length(ttc{ii}{ii}) ;
        % Now label each timepoint in linear indexing
        if ii == 1
            tpid = cat(2, tpid, 1:ntps(ii)) ;
        else
            tpid = cat(2, tpid, addt(ii) + (1:ntps(ii))) ;
        end

        % Store all indices in one array
        tau0extra = fnxtr(tau0(ii), 2) ;
        tauv = ppval(tau0extra, (1:ntps(ii))) ;
        i_tau0j(addt(ii) + (1:ntps(ii)), :) = [ii * ones(ntps(ii), 1), tauv'] ;
    end
    ntp_tot = sum(ntps) ;
    
    % save data
    save(itau0jfn, 'i_tau0j', 'ntp_tot', 'addt', 'ntps') ;
    
    % de-clutter
    clearvars tauv tau0extra tpid ii
end

%% Next build NL, KL, BL of correspondences for master timeline
% NL = neighbor list, KL = spring constant list, and BL = bond list
NLKLBLfn = fullfile(timelineDir, 'NL_KL_BL.mat') ;

if exist(NLKLBLfn, 'file') && ~overwrite
    load(NLKLBLfn, 'NL', 'KL', 'BL') ;
    disp('loaded bond list BL, neighbor list NL, k list KL from disk')
else
    % preallocate Nxlarge array for NL, KL, do not preallocate BL
    NL = zeros(ntp_tot, 30) ;
    KL = zeros(ntp_tot, 30) ;
    first = true ;
    for ii = 1:length(ttc)
        for jj = ii+1:length(ttc)
            % Ensure that there are some correspondences
            if ~isempty(ttc{ii}{jj})
                nodei = round(ttc{ii}{jj}(:, 1)) + addt(ii);
                nodej = round(ttc{ii}{jj}(:, 2)) + addt(jj);
                bladd = [nodei, nodej] ;

                % Add to BL
                if first
                    BL = bladd ;
                    first = false ;
                else
                    BL = cat(1, BL, bladd) ;
                end

                % build NL, KL with redundancy of correspondences built into KL
                for id = 1:length(bladd)
                    pair = bladd(id, :) ;
                    nodei = pair(1) ; 
                    nodej = pair(2) ;
                    % ij
                    if ismember(nodej, NL(nodei, :))
                        ind = find(NL(nodei, :) == nodej) ;
                        assert(NL(nodei, ind) == nodej) ;
                        KL(nodei, ind) = KL(nodei, ind) + 1 ;
                    else
                        firstzero = find(NL(nodei, :)==0, 1) ; 
                        assert(~isempty(firstzero)) 
                        NL(nodei, firstzero) = nodej ; 
                        KL(nodei, firstzero) = 1 ;
                    end

                    % ji
                    if ismember(nodei, NL(nodej, :))
                        ind = find(NL(nodej, :) == nodei) ;
                        assert(NL(nodej, ind) == nodei) ;
                        KL(nodej, ind) = KL(nodej, ind) + 1 ;
                    else
                        firstzero = find(NL(nodej, :)==0, 1) ; 
                        assert(~isempty(firstzero)) 
                        NL(nodej, firstzero) = nodei ; 
                        KL(nodej, firstzero) = 1 ;
                    end
                end
            end
        end
    end
    disp('done building bond list BL, neighbor list NL, k list KL ')
end

%% Build network visualization (pairwise_corr_timeline_XXX.png)

for use_BL = [true false]
    if use_BL
        disp('Visualizing network using BL')
    else
        disp('VIsualizing network using KL, NL')
    end
    for ii = 1:1:length(ttc)    
        close all
        for qq = 1:size(i_tau0j, 1)
            plot(i_tau0j(qq, 2), i_tau0j(qq, 1), '.', 'color', colorset(i_tau0j(qq, 1), :))
            hold on
        end
        % labels
        ylabel('dataset $i$', 'Interpreter', 'Latex')
        xlabel('time, $t_0$', 'Interpreter', 'latex')
        xlim([0, max(i_tau0j(:, 2)) + 1])

        if use_BL
            % Add bonds to plot using BL (non-reciprocal)
            for qq = 1:size(BL, 1)
                if i_tau0j(BL(qq, 1), 1) == ii 
                    plot([i_tau0j(BL(qq, 1), 2), i_tau0j(BL(qq, 2), 2)], ...
                         [i_tau0j(BL(qq, 1), 1), i_tau0j(BL(qq, 2), 1)], '-', ...
                         'color', colorset(i_tau0j(BL(qq, 2), 1), :), 'linewidth', 0.01)
                elseif i_tau0j(BL(qq, 2), 1) == ii
                    plot([i_tau0j(BL(qq, 1), 2), i_tau0j(BL(qq, 2), 2)], ...
                         [i_tau0j(BL(qq, 1), 1), i_tau0j(BL(qq, 2), 1)], '-', ...
                         'color', colorset(i_tau0j(BL(qq, 1), 1), :), 'linewidth', 0.01)
                end
            end

            title(['Time relaxation network: dataset ' num2str(ii)])
            ylim([0, max(i_tau0j(:, 1)) + 1])
            saveas(gcf, fullfile(timelineDir, sprintf('pairwise_corr_timeline_BL_%03d.png', ii)))
            % disp('press any button on figure')
            pause(0.1)
        else
            % Add bonds to plot using NL,KL 
            % cycle only through the linear point indices that belong to ii
            for qq = (addt(ii)+1):(addt(ii) + ntps(ii))
                for nn = NL(qq, :)
                    % plot only if real neighbor
                    if nn > 0
                        plot([i_tau0j(qq, 2), i_tau0j(nn, 2)], ...
                             [i_tau0j(qq, 1), i_tau0j(nn, 1)], '-', ...
                             'color', colorset(i_tau0j(nn, 1), :), 'linewidth', 0.01)
                    end
                end
            end
            title(['Time relaxation network: dataset ' num2str(ii)])
            ylim([0, max(i_tau0j(:, 1)) + 1])
            saveas(gcf, fullfile(timelineDir, sprintf('pairwise_corr_timeline_%03d.png', ii)))
            % disp('press any button')
            pause(0.1)
        end
    end
end
disp('done with network visualization')
close all

%% post-relaxation network visualization (i_tau0j_tau0jrelaxed.mat)

% check for existing result
fn = fullfile(timelineDir, 'i_tau0j_tau0jrelaxed.mat') ;
if exist(fn, 'file') && ~overwrite
    load(fn, 'i_tau0j_tau0jrelaxed')
else
    % Relax: fix the time nodes of the hard reference dataset, move all others
    options = optimset('PlotFcns',@optimplotfval, 'TolFun', 1e-6);
    x0 = i_tau0j(:, 2) ;
    % pop the indices of the fixed times from the array x0
    fixed_ind = find(i_tau0j(:, 1) == hard) ;
    fixed_xx = i_tau0j(fixed_ind, 2) ;
    x0(fixed_ind) = [] ;
    fun = @(x)springEnergy1D(x, BL, fixed_ind, fixed_xx);
    xf = fminsearch(fun, x0, options) ;
    % reinsert indices of fixed times
    xnew1 = xf(1:fixed_ind(1)-1) ;
    xnew2 = xf(fixed_ind(1):end) ;
    xnew = [xnew1; fixed_xx; xnew2 ] ;
    disp('done')

    % Plot the relaxed network result 
    close all
    offset = -0.25 ;
    for qq = 1:size(i_tau0j, 1)
        plot(i_tau0j(qq, 2), i_tau0j(qq, 1) + offset, ...
            'o', 'color', colorset(i_tau0j(qq, 1), :))
        hold on
    end
    for qq = 1:size(i_tau0j, 1)
        plot(xnew(qq), i_tau0j(qq, 1), '.', 'color', colorset(i_tau0j(qq, 1), :))
        hold on
    end

    % Show sliding movement 
    for qq = 1:size(i_tau0j, 1)
        plot([i_tau0j(qq, 2), xnew(qq)], i_tau0j(qq, 1) + [offset, 0], ...
            '-', 'color', colorset(i_tau0j(qq, 1), :))
        hold on
    end

    % labels
    ylabel('dataset $i$', 'Interpreter', 'Latex')
    xlabel('time, $t_0$', 'Interpreter', 'latex')
    xlim([0, max(i_tau0j(:, 2)) + 1])

    title(['Time relaxation network'])
    ylim([0, max(i_tau0j(:, 1)) + 1])
    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [0 0 36 16]) ;
    saveas(gcf, fullfile(outdir, sprintf('relaxation_results.png')))

    % Save the result
    i_tau0j_tau0jrelaxed = cat(2, i_tau0j, xnew) ;
    save(fn, 'i_tau0j_tau0jrelaxed')
end
disp('done')

%% Now look at stripe 7 with standard deviations
% check for existing result
allstripeFn = fullfile(timelineDir, 'allstripes.mat') ;

if exist(allstripeFn, 'file') && ~overwrite
    load(allstripeFn, 'allstripes') ;
else
    % Gather all stripes
    for ii = 1:length(expts)
        disp(['dataset ii = ', num2str(ii)])
        % Load time sequence MIPs of dataset ii
        curvIfn = fullfile(expts{ii}, 'stripe7curves.mat') ;
        disp(['  Loading curves: ' curvIfn])
        load(curvIfn, 'stripe7curves');
        if ii == 1
            allstripes = stripe7curves ;
        else
            allstripes = cat(1, allstripes, stripe7curves) ;
        end
    end
    
    % Save the data allstripes
    save(allstripeFn, 'allstripes')    
end

% Get xlim, ylim for plotting
xmax = 0 ;
xmin = 0 ;
ymin = 0 ;
ymax = 0 ;
for qq = 1:length(allstripes)
    stripe = allstripes{qq} ;
    xmin = min(xmin, min(stripe(:, 1))) ;
    xmax = max(xmax, max(stripe(:, 1))) ;
    ymin = min(ymin, min(stripe(:, 2))) ;
    ymax = max(ymax, max(stripe(:, 2))) ;
end
% Get lateral width
wAP = xmax ;
wDV = ymax ;

%% Back out the mean and variances of leading and trailing curves
statfn = fullfile(corrOutDir, 'curve7stats', 'curve7stats_%04d.mat') ;
datDir = fullfile(corrOutDir, 'curve7stats') ;
stripeTauDir0 = fullfile(corrOutDir, 'aligned_stripes') ;
stripeTauDir = fullfile(stripeTauDir0, 'raw') ;
stripeTauRawDir = fullfile(stripeTauDir, 'raw') ;
stripeTauStatDir = fullfile(stripeTauDir, 'stat') ;
dirs2make = {datDir, stripeTauDir0, stripeTauRawDir, stripeTauStatDir} ;
for di = 1:length(dirs2make)
    dir2make = dirs2make{di} ;
    if ~exist(dir2make, 'dir')
        mkdir(dir2make)
    end
end

% consider each timepoint
window = 0.5 + 1e-5 ;
times2do = ttc{hard}{hard}(:, 1) ;
for qi = 1:length(times2do)
    qq = times2do(qi) ;
    if ~exist(sprintf(statfn, qi), 'file') || overwrite
        close all
        fig = figure('visible', 'off') ;
        disp(['qq = ' num2str(qq)])
        % for this time, get all curves with matching time
        inds = find(abs(i_tau0j_tau0jrelaxed(:, 2) - qq) < window) ;
        stripes = allstripes(inds) ;

        for sn = 1:length(stripes)  
            stripe = stripes{sn} ;        
            % Chop into two curves: leading, trailing
            segs = cutCurveIntoLeadingTrailingSegments(stripe) ;

            % Take mean and variance for leading, trailing
            if sn == 1
                leading = segs{1} ;
                trailing = segs{2} ;
            else
                leading = cat(1, leading, segs{1}) ;
                trailing = cat(1, trailing, segs{2}) ;    
            end
        end
        % swap columns and rows
        lswap = leading(:,[2 1]);
        tswap = trailing(:,[2 1]);
        lstat = binnedstats(lswap) ;    
        % % enforce edges to span the whole image for trailing edge
        % if qi == 1
        %     edges_fixed = (min(lstat(:, 1)) - 0.5): dx : (max(lstat(:, 1))+0.5) ;
        % end
        tstat = binnedstats(tswap) ;

        % Plot the curves
        plot(leading(:, 1)/ wAP, leading(:, 2)/ wDV, '.', 'color', blue)
        hold on;
        plot(trailing(:, 1)/ wAP, trailing(:, 2)/ wDV, '.', 'color', red)
        plot(lstat(:, 2)/ wAP, lstat(:, 1)/ wDV, 'k.')
        plot(tstat(:, 2)/ wAP, tstat(:, 1)/ wDV, 'k.')

        % Save figure
        daspect([wAP, wDV, 1])
        xlim([0 1])
        ylim([0 1])
        title(['stripe 7, t=' num2str(qq)])
        xlabel('AP position [x/L]')
        ylabel('DV position [\phi/C]')
        outfn = fullfile(stripeTauRawDir, sprintf('raw_%04d.png', qq)) ;
        disp(['saving ' outfn])
        saveas(fig, outfn)

        % Convert to fractional values 
        % Get full range of lateral positions (assume contained in first timepoint)
        if qi == 1
            xx = lstat(:, 1) ;
        end
        inds = ismember(xx, lstat(:, 1)) ;
        yy = -ones(size(xx)) ;
        zz = -ones(size(xx)) ;
        yy(inds) = lstat(:, 2) ;
        zz(inds) = lstat(:, 3) ;
        yy(yy<0) = NaN ;
        zz(zz<0) = NaN ;
        leadfrac = cat(2, xx / wDV, yy / wAP, zz / wAP, sqrt(zz) / wAP) ;

        % trailing
        inds = ismember(xx, tstat(:, 1)) ;
        yy = -ones(size(xx)) ;
        zz = -ones(size(xx)) ;
        yy(inds) = tstat(:, 2) ;
        zz(inds) = tstat(:, 3) ;
        yy(yy<0) = NaN ;
        zz(zz<0) = NaN ;
        trailfrac = cat(2, xx / wDV, yy / wAP, zz/ wAP, sqrt(zz) / wAP) ;

        % Save statistics
        leading = lstat ;
        trailing = tstat ; 
        save(sprintf(statfn, qi), 'leading', 'trailing', 'leadfrac', 'trailfrac') ;
        clear leading trailing lstat tstat stripes

        % Check
        close all
        fig = figure('visible', 'off') ;
        plot(leadfrac(:, 2), leadfrac(:, 1), '.'); hold on
        plot(leadfrac(:, 2)- leadfrac(:, 4), leadfrac(:, 1), '.')
        plot(leadfrac(:, 2)+ leadfrac(:, 4), leadfrac(:, 1), '.')
        % Save figure
        daspect([wAP, wDV, 1])
        xlim([0 1])
        ylim([0 1])
        title(['stripe 7, t=' num2str(qq)])
        xlabel('AP position [x/L]')
        ylabel('DV position [\phi/C]')
        outfn = fullfile(stripeTauStatDir, sprintf('frac_%04d.png', qq)) ;
        disp(['saving ' outfn])
        saveas(fig, outfn)
    else
        disp(['tp ' num2str(qq) ' already done'])
    end
end
disp('done with mean and variances of leading and trailing curves')

%% Smooth mean curves in master time
curv7MatFn = fullfile(corrOutDir, 'curve7stats_filtered.mat') ;
if exist(curv7MatFn, 'file')
    disp('Loading leading and trailing stripe7 matrices from disk')
    load(curv7MatFn)
else
    % Grab xx from first timepoint's lstat
    load(sprintf(statfn, 1), 'leading') ;
    xx = leading(:, 1) ;

    % Preallocate matrices of positions and variances in position of leading
    % and trailing stripes
    LX = zeros(length(xx), length(times2do)) ;
    LV = zeros(length(xx), length(times2do)) ;
    LS = zeros(length(xx), length(times2do)) ;
    TX = zeros(length(xx), length(times2do)) ;
    TV = zeros(length(xx), length(times2do)) ;
    TS = zeros(length(xx), length(times2do)) ;
    for qi = 1:length(times2do)
        disp(['Collating time ' num2str(qi)])
        close all
        fig = figure('visible', 'off') ;
        qq = times2do(qi) ;
        load(sprintf(statfn, qi), 'leadfrac', 'trailfrac') ;
        LX(:, qi) = leadfrac(:, 2) ; % position
        LV(:, qi) = leadfrac(:, 3) ; % variance
        LS(:, qi) = leadfrac(:, 4) ; % stdev
        TX(:, qi) = trailfrac(:, 2) ; % position
        TV(:, qi) = trailfrac(:, 3) ; % variance
        TS(:, qi) = trailfrac(:, 4) ; % stdev
    end

    % Now smooth along time dim
    hh = [1 2 3 2 1] / sum([1 2 3 2 1]) ;
    LXs = LX ; 
    LVs = LV ; 
    LSs = LS ;
    TXs = TX ; 
    TVs = TV ; 
    TSs = TS ;
    % Filter each row BEFORE NANs only
    for xi = 1:size(LX, 1) 
        disp(['filtering row ' num2str(xi)])
        % Leading
        xrow = LX(xi, :) ;
        vrow = LV(xi, :) ;
        srow = LS(xi, :) ;
        last = find(isnan(xrow), 1) ;
        rsm = imfilter(xrow(1:last-1), hh, 'replicate') ;
        vsm = imfilter(vrow(1:last-1), hh, 'replicate') ;
        ssm = imfilter(srow(1:last-1), hh, 'replicate') ;
        LXs(xi, 1:last-1) = rsm ;
        LVs(xi, 1:last-1) = vsm ;
        LSs(xi, 1:last-1) = ssm ;
        clear rsm vsm ssm 

        % Trailing
        xrow = TX(xi, :) ;
        vrow = TV(xi, :) ;
        srow = TS(xi, :) ;
        rsm = imfilter(xrow(1:last-1), hh, 'replicate') ;
        vsm = imfilter(vrow(1:last-1), hh, 'replicate') ;
        ssm = imfilter(srow(1:last-1), hh, 'replicate') ;
        TXs(xi, 1:last-1) = rsm ;
        TVs(xi, 1:last-1) = vsm ;
        TSs(xi, 1:last-1) = ssm ;
    end
    % Save to disk
    save(curv7MatFn, 'LX', 'LV', 'LS', 'LXs', 'LVs', 'LSs', 'TX', 'xx')

    % Save the figure versions
    todo = {LX, LS, TX, TS, LXs, LSs, TXs, TSs} ;
    labels = {'AP position of stripe 7 (raw)', ...
        'AP stdev of stripe 7 (raw)', ...
        'AP position of trailing edge stripe 7 (raw)', ...
        'AP stdev of trailing edge stripe 7 (raw)', ...
        'AP position of stripe 7 (smoothed)', ...
        'AP stdev of stripe 7 (smoothed)', ...
        'AP position of trailing edge stripe 7 (smoothed)', ...
        'AP stdev of trailing edge stripe 7 (smoothed)'} ;
    names = {'posL_raw', 'stdL_raw', 'posT_raw', 'stdT_raw', ...
        'posL_smooth', 'stdL_smooth', 'posT_smooth', 'stdT_smooth'} ;
    for i = 1:length(todo)
        close all
        imagesc(todo{i})
        xlabel('master time')
        ylabel('DV position, [\phi/C]')
        title(labels{i})
        cb = colorbar() ;
        if mod(i, 2) == 1
            ylabel(cb, 'AP position') 
        else
            ylabel(cb, '\sigma_{x/L}')
        end
        saveas(gcf, fullfile(corrOutDir, ['position_heatmap' names{i} '.png'])) ;
    end
end
disp('done with time filtering')

%% Minimize for lateral motion now that we have averages in master time
% For each master time, rotate each curve to match average
window = 0.5 + 1e-5 ;
times2do = ttc{hard}{hard}(:, 1) ;
datDir = fullfile(corrOutDir, 'curve7stats_collapsed') ;
statCollapseFn = fullfile(datDir, 'curve7stats_collapsed_%04d.mat') ;
stripeTauCollDir = fullfile(corrOutDir, 'aligned_stripes', 'collapsed') ;
stripeTauCollRawDir = fullfile(stripeTauCollDir, 'raw') ;
stripeTauCollStatDir = fullfile(stripeTauCollDir, 'stat') ;
dirs2make = {datDir, stripeTauCollDir, ...
    stripeTauCollRawDir, stripeTauCollStatDir} ;
for di = 1:length(dirs2make)
    dir2make = dirs2make{di} ;
    if ~exist(dir2make, 'dir')
        mkdir(dir2make)
    end
end

% Check if this has already been done
redo_optimization = false ;
for qi = 1:length(times2do)
    datfn = sprintf(statCollapseFn, qi) ;
    if ~exist(datfn, 'file')
        redo_optimization = true ;
    end
end

% Consider all times
if redo_optimization || overwrite 
    for qi = 1:length(times2do)
        close all
        qq = times2do(qi) ;
        disp(['qq = ' num2str(qq)])
        % for this time, get all curves with matching time
        inds = find(abs(i_tau0j_tau0jrelaxed(:, 2) - qq) < window) ;
        stripes = allstripes(inds) ;

        for sn = 1:length(stripes)  
            stripe = stripes{sn} ;        
            % Chop into two curves: leading, trailing
            segs = cutCurveIntoLeadingTrailingSegments(stripe) ;
            lead = segs{1} ;
            trail = segs{2} ;
            lead(:, 1) = lead(:, 1) / wAP ;
            lead(:, 2) = lead(:, 2) / wDV ;
            trail(:, 1) = trail(:, 1) / wAP ;
            trail(:, 2) = trail(:, 2) / wDV ;
            % swap xy
            lead(:, [1 2]) = lead(:, [2 1]) ;
            trail(:, [1 2]) = trail(:, [2 1]) ;

            % Optimize for DV motion, add to list of leading for this TP
            x0 = [0., 0.] ;
            refcurv = zeros(length(LXs(:, qi)), 2) ;
            refcurv(:, 1) = double(1:length(LXs(:, qi))) / double(length(LXs(:, qi))) ;
            refcurv(:, 2) = LXs(:, qi) ;
            fun = @(x)ssrCurves(lead + [x(1) x(2)], refcurv, true, true) ;
            shifts = fminsearch(fun, x0) ;
            leado = lead + [shifts(1), shifts(2)] ;
            % Modulo for wDV
            leado(:, 1) = mod(leado(:, 1), 1) ;

            % % check it
            % plot(lead(:, 1), lead(:, 2), '.'); hold on;
            % plot(leado(:, 1), leado(:, 2), '.'); hold on;
            % plot(refcurv(:, 1), refcurv(:, 2), '.')

            % trailing optimization
            refcurv = zeros(length(TXs(:, qi)), 2) ;
            refcurv(:, 1) = double(1:length(TXs(:, qi))) / double(length(TXs(:, qi))) ;
            refcurv(:, 2) = TXs(:, qi) ;
            fun = @(x)ssrCurves(trail + [x(1) x(2)], refcurv, true, true) ;
            shifts = fminsearch(fun, x0) ;
            trailo = trail + [shifts(1), shifts(2)] ;
            % Modulo for wDV
            trailo(:, 1) = mod(trailo(:, 1), 1) ;

            % % check it
            % close all
            % plot(trail(:, 1), trail(:, 2), '.'); hold on;
            % plot(trailo(:, 1), trailo(:, 2), '.'); hold on;
            % plot(refcurv(:, 1), refcurv(:, 2), '.')
            % pause(1)

            if sn == 1
                leading = leado ;
                trailing = trailo ;
            else
                leading = cat(1, leading, leado) ;
                trailing = cat(1, trailing, trailo) ;
            end
        end
        % Do stats
        dx = 1.0 / size(xx, 1) ;
        edges = 0:dx:1.0 ;
        lstat = binnedstats(leading, edges) ;
        tstat = binnedstats(trailing, edges) ;

        % Plot the curves
        close all
        fig = figure('visible', 'off') ;
        plot(leading(:, 2), leading(:, 1), '.', 'color', blue)
        hold on;
        plot(trailing(:, 2), trailing(:, 1), '.', 'color', red)
        plot(lstat(:, 2), lstat(:, 1), 'k.')
        plot(tstat(:, 2), tstat(:, 1), 'k.')

        % Save figure
        daspect([wAP, wDV, 1])
        xlim([0 1])
        ylim([0 1])
        title(['collapsed stripe 7, t=' num2str(qq)])
        xlabel('AP position [x/L]')
        ylabel('DV position [\phi/C]')
        outfn = fullfile(stripeTauCollRawDir, ...
            sprintf('collapsed_raw_statcurve7_%04d.png', qq)) ;
        disp(['Saving ' outfn])
        saveas(fig, outfn)

        % Convert to matrix-style curves (all same width, same sampling)
        % Get full range of lateral positions (assume contained in first timepoint)
        if qi == 1
            xx = lstat(:, 1) ;
        end

        % match indices
        % [~, inds] = min(abs(pdist2(xx, lstat(:, 1))), [], 2) ;
        inds = false(size(xx)) ;
        for dmyk = 1:length(inds)
            if any(abs(xx(dmyk) - lstat(:, 1)) < dx * 0.6)
                inds(dmyk) = true ;
            end
        end
        yy = -ones(size(xx)) ;
        zz = -ones(size(xx)) ;
        yy(inds) = lstat(:, 2) ;
        zz(inds) = lstat(:, 3) ;
        yy(yy<0) = NaN ;
        zz(zz<0) = NaN ;
        leadfrac = cat(2, xx, yy, zz, sqrt(zz)) ;

        % trailing
        % [~, inds] = min(abs(pdist2(xx, tstat(:, 1))), [], 2) ;
        inds = false(size(xx)) ;
        for dmyk = 1:length(inds)
            if any(abs(xx(dmyk) - tstat(:, 1)) < dx * 0.5)
                inds(dmyk) = true ;
            end
        end
        yy = -ones(size(xx)) ;
        zz = -ones(size(xx)) ;
        yy(inds) = tstat(:, 2) ;
        zz(inds) = tstat(:, 3) ;
        yy(yy<0) = NaN ;
        zz(zz<0) = NaN ;
        trailfrac = cat(2, xx, yy, zz, sqrt(zz)) ;

        % Save statistics
        leading = lstat ;
        trailing = tstat ; 
        datfn = sprintf(statCollapseFn, qi) ;
        disp(['Saving data to ' datfn])
        save(datfn, 'leading', 'trailing', 'leadfrac', 'trailfrac') ;
        clear leading trailing lstat tstat stripes

        % Check
        close all
        fig = figure('visible', 'off') ;
        plot(leadfrac(:, 2), leadfrac(:, 1), '.'); hold on
        plot(leadfrac(:, 2) - leadfrac(:, 4), leadfrac(:, 1), '.')
        plot(leadfrac(:, 2) + leadfrac(:, 4), leadfrac(:, 1), '.')
        % Save figure
        daspect([wAP, wDV, 1])
        xlim([0 1])
        ylim([0 1])
        title(['collapsed stripe 7, t=' num2str(qq)])
        xlabel('AP position [x/L]')
        ylabel('DV position [\phi/C]')
        outfn = fullfile(stripeTauCollStatDir, ...
            sprintf('collapsed_statcurve7_%04d.png', qq)) ;
        disp(['Saving ' outfn])
        saveas(fig, outfn)
    end
    disp('Done with collapsing curves')
else
    disp('Already did collapse of curves')
end
disp('done with curve collapse')

%% Smooth collapsed mean curves in master time
filtWidth = 10;
filtSigma = 6 ;
imageFilter = fspecial('gaussian',filtWidth,filtSigma);
        
c7MatFiltFn = fullfile(corrOutDir, 'curve7stats_collapsed_filtered.mat') ;
if ~exist(c7MatFiltFn, 'file') || overwrite 
    % Grab xx from first timepoint's lstat 
    load(sprintf(statCollapseFn, 1), 'leading') ;
    xx = leading(:, 1) ;

    % Preallocate matrices of positions and variances in position of leading
    % and trailing stripes
    LX = nan(length(xx), length(times2do)) ;
    LV = nan(length(xx), length(times2do)) ;
    LS = nan(length(xx), length(times2do)) ;
    TX = nan(length(xx), length(times2do)) ;
    TV = nan(length(xx), length(times2do)) ;
    TS = nan(length(xx), length(times2do)) ;
    for qi = 1:length(times2do)
        disp(['Collating time ' num2str(qi)])
        close all
        fig = figure('visible', 'off') ;
        qq = times2do(qi) ;
        load(sprintf(statCollapseFn, qi), 'leadfrac', 'trailfrac') ;
        LX(:, qi) = leadfrac(:, 2) ; % position
        LV(:, qi) = leadfrac(:, 3) ; % variance
        LS(:, qi) = leadfrac(:, 4) ; % stdev
        TX(:, qi) = trailfrac(:, 2) ; % position
        TV(:, qi) = trailfrac(:, 3) ; % variance
        TS(:, qi) = trailfrac(:, 4) ; % stdev
    end
    LX(LX==0) = NaN ;
    LV(LV==0) = NaN ;
    LS(LS==0) = NaN ;
    TX(TX==0) = NaN ;
    TV(TV==0) = NaN ;
    TS(TS==0) = NaN ;
    LX = (fillmissing(LX', 'previous'))' ;
    LV = (fillmissing(LV', 'previous'))' ;
    LS = (fillmissing(LS', 'previous'))' ;
    TX = (fillmissing(TX', 'previous'))' ;
    TV = (fillmissing(TV', 'previous'))' ;
    TS = (fillmissing(TS', 'previous'))' ;

    % Now smooth along time dim
    hh = [1 2 3 4 5 6 5 4 3 2 1] / sum([1 2 3 4 5 6 5 4 3 2 1]) ;
    LXs = LX ; 
    LVs = LV ; 
    LSs = LS ;
    TXs = TX ; 
    TVs = TV ; 
    TSs = TS ;
    % Filter each row BEFORE NANs only, take sliding max with nans set to zero
    for xi = 1:size(LX, 1) 
        disp(['filtering row ' num2str(xi)])
        % Leading
        xrow = LX(xi, :) ;
        vrow = LV(xi, :) ;
        srow = LS(xi, :) ;
        last = find(isnan(xrow), 1) ;
        rsm = imfilter(xrow(1:last-1), hh, 'replicate') ;
        vsm = slidingmax(vrow, length(hh)) ;
        vsm = imfilter(vsm, hh, 'replicate') ;
        ssm = slidingmax(srow, length(hh)) ;
        ssm = imfilter(ssm, hh, 'replicate') ;
        LXs(xi, 1:last-1) = rsm ;
        LVs(xi, :) = vsm ;
        LSs(xi, :) = ssm ;
        clear rsm vsm ssm 

        % Trailing
        xrow = TX(xi, :) ;
        vrow = TV(xi, :) ;
        srow = TS(xi, :) ;
        rsm = imfilter(xrow(1:last-1), hh, 'replicate') ;
        vsm = slidingmax(vrow, length(hh), true) ;
        vsm = imfilter(vsm, hh, 'replicate') ;
        ssm = slidingmax(srow, length(hh), true) ;
        ssm = imfilter(ssm, hh, 'replicate') ;
        TXs(xi, 1:last-1) = rsm ;
        TVs(xi, :) = vsm ;
        TSs(xi, :) = ssm ;
    end
    
    % Smooth in 2D twice
    LVs = nanconv(LVs, imageFilter, 'nanout');
    LSs = nanconv(LSs, imageFilter, 'nanout');
    TVs = nanconv(TVs, imageFilter, 'nanout');
    TSs = nanconv(TSs, imageFilter, 'nanout');
    
    % Save to disk
    save(fullfile(corrOutDir, 'curve7stats_collapsed_filtered.mat'), ...
        'LX', 'LV', 'LS', 'LXs', 'LVs', 'LSs', ...
        'TX', 'TV', 'TS', 'TXs', 'TVs', 'TSs', 'xx')

    % Save the figure versions
    todo = {LX, LS, TX, TS, LXs, LSs, TXs, TSs} ;
    labels = {'AP position of stripe 7 (raw)', ...
        'AP stdev of stripe 7 (raw)', ...
        'AP position of trailing edge stripe 7 (raw)', ...
        'AP stdev of trailing edge stripe 7 (raw)', ...
        'AP position of stripe 7 (smoothed)', ...
        'AP stdev of stripe 7 (smoothed)', ...
        'AP position of trailing edge stripe 7 (smoothed)', ...
        'AP stdev of trailing edge stripe 7 (smoothed)'} ;
    names = {'posL_raw', 'stdL_raw', 'posT_raw', 'stdT_raw', ...
        'posL_smooth', 'stdL_smooth', 'posT_smooth', 'stdT_smooth'} ;
    for i = 1:length(todo)
        close all
        imagesc(todo{i})
        xlabel('master time')
        ylabel('DV position, [\phi/C]')
        title(labels{i})
        cb = colorbar() ;
        if mod(i, 2) == 1
            ylabel(cb, 'AP position') 
        else
            ylabel(cb, '\sigma_{x/L}')
        end
        saveas(gcf, fullfile(corrOutDir, ['position_heatmap_collapsed_' names{i} '.png'])) ;
    end
else
    disp('Filtered curve7 matrix already on disk')
end
disp('done with final collapsed output')

%% How accurate can we measure the timing?
calibDir = fullfile(alignDir, 'alignment_calibration') ;
load(fullfile(corrOutDir, 'curve7stats_collapsed_filtered.mat'), ...
    'LXs', 'LVs', 'xx') ;

hands_on = true ;
smooth_var = false ;
save_snaps = false ;
refcurvX = xx ;
refcurvsY = LXs' ;
refvars = LVs' ;
% Align Stripe 7 for each frame with chisquare analysis
close all
fig = figure('visible', 'off') ;
for exptii = 4:length(expts)
    clf
    disp(['dataset ii = ', num2str(exptii)])
    disp(expts{exptii})
    
    % create dir to put snapshots
    figdir = fullfile(calibDir, sprintf('dataset_%04d_snaps', exptii)) ;
    if ~exist(figdir, 'dir')
        mkdir(figdir)
    end

    % Load time sequence MIPs of dataset ii
    curvIfn = fullfile(expts{exptii}, 'stripe7curves.mat') ;
    load(curvIfn, 'stripe7curves') ;
    chisq = zeros([length(stripe7curves), 1]) ;
    
    % preallocate
    tmatches = zeros(length(stripe7curves), 1) ;
    tuncs = zeros(length(stripe7curves), 1) ;
    tmatches_ssr = zeros(length(stripe7curves), 1) ;
    tuncs_ssr = zeros(length(stripe7curves), 1) ;
    tt = 1:length(stripe7curves) ;
    for qq = 1:10:length(tt)
        if mod(qq, 10) == 0
            disp(['calibration: qq = ' num2str(qq)])
        end
        stripe = stripe7curves{qq} ;
        xx = stripe(:, 1) / wAP ;
        yy = stripe(:, 2) / wDV ;
        segs = cutCurveIntoLeadingTrailingSegments([xx, yy]) ;
        lead = segs{1} ;
        curv = lead(:, [2 1]) ;
        [chisq, chisqn, ssr] = chisquareCurves(curv, refcurvX, ...
            refcurvsY, refvars, smooth_var, true) ;
        % error('here')
        chisqn = fillmissing(chisqn, 'movmedian', 5) ;
        if hands_on
            [tmatch, unc, fit_coefs, nidx, pidx] = ...
                chisqMinUncInteractiveDomain(chisqn, 4, 10, ssr) ;
        else        
            [tmatch, unc, fit_coefs] = chisqMinUncertainty(chisqn, 4, 10) ;
            nidx = round(tmatch) - 10 ;
            pidx = round(tmatch) + 10 ;    
        end
        % Plot the result
        if save_snaps
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
            xlabel(['time [min], dataset' num2str(exptii)])
            title(['dataset ' num2str(exptii) ': $t=$' sprintf('%04d', qq)], ...
                'Interpreter', 'Latex')
            figfn = fullfile(figdir, sprintf('chisq_snap_%04d.png', qq)) ;
            disp(['Saving ' figfn])
            saveas(fig, figfn)
            clf
        end
        
        % replace NaNs with large dt
        unc(isnan(unc)) = 50 ;
        tmatches(qq) = tmatch ;
        tuncs(qq) = unc ;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SSR 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Filter to get minimum of SSR
        % windowWidth = 3;  % for quick smoothing
        % kernel = ones(windowWidth, 1) / windowWidth;
        % ssdsm = filter(kernel, 1, ssds);
        nfit = pidx - nidx ;
        [matchtime, matchtime_unc] = matchTimeSSR(ssr, nfit) ;
        
        % store ssr timematches
        tmatches_ssr(qq) = matchtime ;
        tuncs_ssr(qq) = matchtime_unc ;
        
        
        close all
    end
    
    % Plot best fit with errors on cij plot with hard reference
    ijstr = [ '_%02d_%02d' extn ] ;
    cfn = fullfile(corrDatOutDir, sprintf(['corr' ijstr '.mat'], exptii, hard)) ;
    load(cfn, 'cij')
    close all
    fig = figure('visible', 'off') ;
    imagesc(cij)
    hold on;
    plot(tmatches, tt, 'o-', 'color', 'k')
    plot(tmatches - tuncs, tt, '--', 'color', 'k')
    plot(tmatches + tuncs, tt, '--', 'color', 'k')
    plot(tmatches_ssr, tt, 'o', 'color', 'r')
    ylabel(['time, dataset ' num2str(exptii) ])
    xlabel(['time, dataset ' num2str(hard) ])
    cb = colorbar() ;
    ylabel(cb, 'cross correlation')
    figoutfn = fullfile(calibDir, ...
        sprintf('dset%02d_calibration_cij_tmatches.png', exptii));
    disp(['Saving ' figoutfn])
    saveas(fig, figoutfn)
    % with SSR
    plot(tmatches_ssr - tuncs_ssr, tt, '--', 'color', 'r')
    plot(tmatches_ssr + tuncs_ssr, tt, '--', 'color', 'r')
    figoutfn = fullfile(calibDir, ...
        sprintf('dset%02d_calibration_cij_tmatches_withSSR.png', exptii));
    disp(['Saving ' figoutfn])
    saveas(fig, figoutfn)
    
    % View it
    set(gcf, 'visible', 'on')
    pause(1)
    close all
    
end


