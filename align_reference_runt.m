% align_reference_runt.m 
% Script for aligning dynamic runt nanobody data against each other
% 
% NPMitchell 2020

%% ADD PATHS
basicsDir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/gut_matlab/basics/' ;
tiffDir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/gut_matlab/tiff_handling/' ;
plottingDir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/gut_matlab/plotting/' ;
codeDir = '/Users/npmitchell/Box/Flies/code/time_alignment_2020/' ;
addpath(basicsDir) ;
addpath(tiffDir) ;
addpath(plottingDir) ;
addpath(codeDir) ;

%% OPTIONS
% Save each runt nanobody (curated) MIP as ./date/cylinder1_max.tif
runtNBodyDir = './Runt-Nanobody/' ;
mipfn = 'cylinder1_max.tif' ;
outdir = './alignment' ;
ssfactor = 4 ;              % subsampling factor before computing corr
% Correlation options
corr_method = 'phase' ; % realspace or phase method for correlation 
                            % If realspace, does not tranlate, phase allows
                            % dx,dy then computes realspace corr on shifted
                            % image.
stripe7corr_method

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
                        'Backspace=manual ROI, a=anteriorCC, e=erode, d=dilate, r=reuse prev ROI'] ;
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
                        % Clear and reset to autogen
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
        cfn = fullfile(outdir, sprintf(['corr' ijstr '.mat'], ii, jj)) ;
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
                        sz1 = min(size(imI, 1), size(imJ, 1)) ;
                        sz2 = min(size(imI, 2), size(imJ, 2)) ;
                        
                        IR = zeros(sz1, sz2) ;
                        JR = zeros(sz1, sz2) ;
                        IR() ;
                        alignImageX()
                        cij(ti, tj) = corr2(imI, imJ) ;
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
            cfn = fullfile(outdir, sprintf(['cij' ijstr '.png'], ii, jj)) ;
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
                cfn = fullfile(outdir, sprintf(['dxij' ijstr '.png'], ii, jj)) ;
                cb = colorbar() ;
                ylabel(cb, 'correlation')
                set(gca,'YDir','normal')
                axis tight
                title('phase correlation offset \deltax')
                saveas(gcf, cfn)        
                % delta y 
                clf
                imagesc(dxij)
                xlabel(['timepoint for dataset ' num2str(jj)])
                ylabel(['timepoint for dataset ' num2str(ii)])
                axis equal
                cfn = fullfile(outdir, sprintf(['dyij' ijstr '.png'], ii, jj)) ;
                cb = colorbar() ;
                ylabel(cb, 'correlation')
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
                        d = pdist2(cI, cJ);
                        [mindist, idx] = min(d, [], 2);
                        cij(ti, tj) = mean(mindist) ;
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
        else
            disp('Skipping since already done')
        end
        
    end
end


%% Minimize difference between timing points
% {tj} = {dtj} + frame, where 
% Fit line of best fit to each ridge extraction

