%% Examine gradients in each label for all dirs in WT
% NPMitchell 2019
% todo: write determineRuntTime() function to map Runt pattern onto
% morphological time axis.
%
% genoDirs are WT, Twist, Ftz, etc
% genoDirs are things like "Eve_Runt_Nero" within the genoDir
% Each genoDir contains many exptDirs labeled by, for ex, the date.
% Each exptDir contains pullbacks with label indices in same order as
%   the genoDir name ("Eve_Runt_Neuro" must contain 3 label pullbacks).
% 
close all
clear

%% Options
view = true ;              % view intermediate results
save_ims = true ;           % save images of the gradients
overwrite_images = false ;   % overwrite the gradient and smoothed images
magfactor = 10 ;            % Factor to multiply gradients to uniformly increase brightness
gather_timebins = false ;    % Average images within time bins
overwrite_binims = false ;  % overwrite the gradient and smoothed images for gathered timebins
rootDir = '../' ;           
genoDirs = {'WT'} ;       % Which directories (genotypes) to analyze
prepend = 'Max_Cyl*_2_000001_c' ;       % string before label index
exten = '.tif' ;    % string after label index
timematfn = 'timematch_EveRunt_tmin27_tmax45.mat' ;
timematfn2 = 'timematch_eve_tmin27_tmax45.mat' ;
% Declare what scales to examine (in pixels)
sigmas = [5, 10, 20, 30, 50] ;
steps = [1] ;
selectedLabels = {'Runt', 'Slp', 'Ftz', 'Paired', 'Eve'} ;
cdf_min = 0.01 ;
cdf_max = 0.999 ;

%% Examine each mutation directory, starting with WT
for kk = 1:length(genoDirs)
    
    % Establish I/O dirs
    genoDir = fullfile(rootDir, genoDirs{kk}) ;
    outDatDir = [genoDir filesep] ;
    
    % Generate folders for output
    if ~exist(outDatDir, 'dir')
        mkdir(outDatDir)
    end

    %% Build the lookuptable
    a = lookupTable ;
    a = a.buildLookup('') ;
    
    %% For each labelDir, look at each embryo
    for chii = 1:length(allLabels)
        label = allLabels{chii} ;
        % Make an output directory for this label only in the genoDir
        outFigDir = fullfile(genoDir, ['figures' filesep label filesep]) ;

        if ~exist(outFigDir, 'dir')
            mkdir(outFigDir)
        end
        
        % Select selected Channels: debug
        if any(contains(selectedLabels, label))
            lut = a.map(label) ;
            disp('done building fileNames for this label')
    
            fileNames = lut.fileNames ;
            
            %% Now load all filenames and make hdf5 dataset 
            % Examine on all scales in kernel and gradient
            for ii = 1:length(sigmas)
                sigma = sigmas(ii) ;
                disp(['Smoothing with sigma = ' num2str(sigma)])

                %% Smooth the data with a gaussian filter and store everything in a struct
                % The struct "data" contains:
                %       im : the raw images 
                %       lim : the intensity limits
                %       scale : the rescaled data according to the intensity limits
                %       smooth : the filtered data images

                for iii = 1:length(fileNames)
                    disp(fileNames(iii))
                end
                
                nFiles = length(fileNames);
                kernel = 5*sigma;
                data = struct('im',[],'lim',[],'scale',[],'smooth',[]);

                %% Load all files
                % fbar = waitbar(0, '') ;
                for k = 1:nFiles
                    msg = ['Reading file ' num2str(k) '/' num2str(nFiles)];
                    % waitbar(k/nFiles, fbar, msg)
                    im = imread(fileNames{k});
                    [f,x] = ecdf(im(:));
                    f1 = find(f>cdf_min, 1, 'first');
                    f2 = find(f<cdf_max, 1, 'last');
                    %imshow(tmp,[x(f1) x(f2)])
                    lim = [x(f1) x(f2)];
                    scale = mat2gray(double(im),double(lim));
                    %figure,imshow(scale)
                    smoothim = imfilter(scale, fspecial('gaussian',kernel,sigma));
                    
                    %% Iterate over all gradient steps
                    for jj = 1:length(steps)
                        step = steps(jj) ;
                        disp(['im ' num2str(k) ' > Gradient with step ' num2str(step)])
                        % Make string representing the gradient step and the smoothing sigma
                        % extensm = sprintf('_sigma%03d', sigma) ;
                        sigstepstr = sprintf('sigma%03d_step%03d', sigma, step) ;
                        exten = ['_' sigstepstr] ;
                        
                        % Take gradients of the average of all data images
                        outFigStepDir = fullfile(fileDirs{k}, sigstepstr) ;
                        sdir = fullfile(outFigStepDir, 'smooth') ;
                        mdir = fullfile(outFigStepDir, 'gradient_magnitude') ;
                        xdir = fullfile(outFigStepDir, 'gradient_dx_magnitude') ;
                        ydir = fullfile(outFigStepDir, 'gradient_dy_magnitude') ;
                        
                        dirs2make = {sdir, mdir, xdir, ydir} ;
                        for qqq=1:length(dirs2make)
                            dir2make = dirs2make{qqq};
                            if ~exist(dir2make, 'dir')
                                mkdir(dir2make)
                            end
                        end
                        % Name the output images
                        sfn = fullfile(sdir, [label '_smooth' exten '.tif']) ;
                        mfn = fullfile(mdir, [label '_gradient_magnitude' exten '.tif']) ;
                        xfn = fullfile(xdir, [label '_gradient_dx_magnitude' exten '.tif']) ;
                        yfn = fullfile(ydir, [label '_gradient_dy_magnitude' exten '.tif']) ;
                        sexist = exist(sfn, 'file') ;
                        mexist = exist(mfn, 'file') ;
                        xexist = exist(xfn, 'file') ;
                        yexist = exist(yfn, 'file') ;
                        imsexist = sexist && mexist && xexist && yexist ;
                        if save_ims && (~imsexist || overwrite_images)
                            % Take gradients of the average of all data images
                            [X,Y] = meshgrid(1:size(smoothim, 1), 1:size(smoothim, 2));
                            %surf(X,Y,double(imavg)'), shading interp
                            [dX,dY] = gradient(double(smoothim(1:end, 1:end))', step);
                            dX = magfactor * dX ;
                            dY = magfactor * dY ;

                            gradmag = sqrt(dX.^2+dY.^2);
                            if view
                                close all
                                figure,
                                %quiver(X(1:step:end,1:step:end),Y(1:step:end,1:step:end),dX,dY)
                                pcolor(Y,X,gradmag), shading interp, axis equal, %colormap gray
                                title([label ' gradient, $|\nabla I|$ '], 'Interpreter', 'Latex')

                                figure, 
                                pcolor(Y,X,abs(dX)), shading interp, axis equal
                                title([label ' $x$ gradient, $|\partial_x I|$ '], 'Interpreter', 'Latex')

                                figure, 
                                pcolor(Y,X,abs(dY)), shading interp, axis equal
                                title([label ' $y$ gradient, $|\partial_y I|$ '], 'Interpreter', 'Latex')
                            end

                            % Save images of the gradients
                            imwrite(smoothim, sfn)
                            imwrite(gradmag', mfn) 
                            imwrite(abs(dX)', xfn)
                            imwrite(abs(dY)', yfn)
                            disp(['saved images for ' label ' ' num2str(k) '/' num2str(nFiles)])
                        end  
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Save image of the composite and variance
                if gather_timebins
                    
                    % Cast each experiment into a timebin
                    edges = linspace(min(fileTimes) - eps, max(fileTimes) + eps, 8) ;
                    timebins = discretize(fileTimes, edges) ;
                    
                    % Consider each bin
                    for qq = 1:max(timebins)
                        % Get all experiments in current bin qq
                        inbin = find(timebins == qq) ;
                        [time, inds] = sort(fileTimes(inbin)) ;
                        time = time / 60 ;  % convert to minutes
                        inbin = inbin(inds) ;
                        
                        % Build smooths compilation
                        rr = 1 ;
                        for pp = inbin
                            disp(['Adding ' fileNames{pp} ' to gathered bin']) 
                            im = imread(fileNames{pp}) ;
                            % fix limits for each layer
                            [f,x] = ecdf(im(:));
                            f1 = find(f>cdf_min, 1, 'first');
                            f2 = find(f<cdf_max, 1, 'last');
                            lim = [x(f1) x(f2)];
                            im = mat2gray(double(im),double(lim));
                            
                            if rr == 1 
                                % Discern the size of the image and create
                                % stack of gathered ims
                                pxsz = size(im, 1) ;
                                pysz = size(im, 2) ;
                                smooths = zeros(length(inbin), pxsz, pysz) ;
                                smooths(rr, :, :) = im ;
                            else
                                % todo: allow slight translation here?
                                
                                % Here we check that images are same size
                                % and add to the stack of gathered ims
                                if all(size(im) == [pxsz, pysz])
                                    smooths(rr, :, :) = im;
                                else
                                    deltax = round(0.5 * size(im, 1) - pxsz) ;
                                    deltay = round(0.5 * size(im, 2) - pysz) ;
                                    if deltax > -0.5
                                        smooths(rr, :, :) = im((1+deltax):(pxsz+deltax), (1+deltay):(pysz+deltay)) ;
                                    else
                                        smooths(rr, (1+deltax):(pxsz-deltax), (1+deltay):(pysz-deltay)) = im ;
                                    end
                                end
                            end
                            rr = rr + 1 ;
                        end
                        
                        % Check it
                        % disp('Use the right and left arrows to scan through the layers')
                        % disp('--> Press Enter when the first (lowest index) layer is found')
                        % flipThroughStackFindLayer(smooths, 'Arrows to scan, <Enter> for min', 1) ; 

                        % Gather smoothed images for each time bin 
                        imsmooth_bin = squeeze(mean(smooths, 1)) ;
                        % Fix limits
                        [f,x] = ecdf(imsmooth_bin(:));
                        f1 = find(f>cdf_min, 1, 'first');
                        f2 = find(f<cdf_max, 1, 'last');
                        lim = [x(f1) x(f2)];
                        imsmooth_bin = mat2gray(double(imsmooth_bin),double(lim));
                        
                        % Define extension for filename
                        exten = sprintf('_sigma%03d_timebin%02d', sigma, qq) ;
                        outFigSigDir = fullfile(outFigDir, [sprintf('sigma%04d', sigma) filesep]) ;

                        if ~exist(outFigSigDir, 'dir')
                            mkdir(outFigSigDir)
                        end
                        if save_ims
                            % save the average
                            imwrite(imsmooth_bin, fullfile(outFigSigDir, ['composite' exten '.png']))
                            % save the variance (fractional)
                            imstd = squeeze(std(smooths, 1)) ;
                            fn = fullfile(outFigSigDir, ['composite_std' exten '.png']) ;
                            imwrite(imstd, fn)
                             
                            % % Get all Runt images
                            % if strcmp(label, 'Runt')
                            %     for ij = 1:size(imsmooth, 3)
                            %         fn = fullfile(outFigSigDir, ['ex' sprintf('%03d', ij) '.png']) ;
                            %         imwrite(imsmooth(:, :, ij), fn) ;
                            %     end
                            % end
                        end

                        %% Iterate over all gradient steps
                        for jj = 1:length(steps)
                            step = steps(jj) ;
                            disp(['  > Timebin ' qq ' gradient with step ' num2str(step)])

                            %% Take gradients of the average of all data images
                            [X,Y] = meshgrid(1:size(imsmooth_bin, 1), 1:size(imsmooth_bin, 2));
                            %surf(X,Y,double(imavg)'), shading interp
                            [dX,dY] = gradient(double(imsmooth_bin(1:end, 1:end))', step) ;
                            dX = magfactor * dX ;
                            dY = magfactor * dY ;
                            
                            % Make string representing the gradient step and the smoothing sigma
                            exten = sprintf('_sigma%03d_step%03d_timebin%02d', sigma, step, qq) ;

                            gradmag = sqrt(dX.^2+dY.^2);
                            %rnt_mag = double(imavg)';
                            if view
                                figure,
                                %quiver(X(1:step:end,1:step:end),Y(1:step:end,1:step:end),dX,dY)
                                pcolor(Y,X,gradmag), shading interp, axis equal, %colormap gray
                                title('Runt gradient, $|\nabla I|$ ', 'Interpreter', 'Latex')

                                figure, 
                                pcolor(Y,X,abs(dX)), shading interp, axis equal
                                title('Runt $x$ gradient, $|\partial_x I|$ ', 'Interpreter', 'Latex')

                                figure, 
                                pcolor(Y,X,abs(dY)), shading interp, axis equal
                                title('Runt $y$ gradient, $|\partial_y I|$ ', 'Interpreter', 'Latex')
                            end

                            % Save images of the gradients
                            mfn = fullfile(outFigSigDir, ['gradient_magnitude' exten '.tif']) ;
                            xfn = fullfile(outFigSigDir, ['gradient_dx_magnitude' exten '.tif']) ;
                            yfn = fullfile(outFigSigDir, ['gradient_dy_magnitude' exten '.tif']) ;
                            mexist = exist(mfn, 'file') ;
                            xexist = exist(xfn, 'file') ;
                            yexist = exist(yfn, 'file') ;
                            if save_ims || overwrite_binims
                                imwrite(gradmag', mfn)
                                imwrite(abs(dX)', xfn)
                                imwrite(abs(dY)', yfn)
                            end
                            disp('Done saving gradient images for composite')
                        end
                    end
                    disp('Done gathering timebins')
                end
            end
            disp(['Done with this sigma: ' sigma])
        end
    end
    disp(['Done with this genoDir: ' genoDir])
end
disp('Done with all genoDirs')

        
%% Take gradients of the unsmoothed data and view it
% im_raw = cat(3,data.im);
% rnt_rav = mean(im_raw,3);
% % Normalize it
% rnt_rav_norm = ((rnt_rav-min(min(rnt_rav)))./(max(max(rnt_rav))-(min(min(rnt_rav)))));
% rnt_rav = rnt_rav_norm;
% 
% cropped_rnt = rnt_rav_norm(:,800:1300);
% mean_cropped_rnt = nanmean(cropped_rnt,2);
% plot(mean_cropped_rnt,(1:2050),'s','LineStyle','none')
