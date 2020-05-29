function makeGradientImages(da, selectedLabels, sigmas, steps,...
                            cdf_minmax, overwrite)
    % 
    % Parameters
    % ----------
    % da : dynamicAtlas class instance 
    %   dynamic atlas class for all genotypes
    % selectedLabels : cell of strings or empty (optional, default=[])
    %   the labels to consider, if nonempty
    % sigmas : int array (optional, default=[5, 10, 20, 30, 50])
    %   the widths for smoothing to be applied on gradients and intensity
    %   data
    % steps : int array (optional, default=[1])
    %
    %   
    % cdf_minmax : 
    % overwrite : 
    %
    % Returns
    % -------
    % 
    %
    % NPMitchell 2020
    
    % Declare what scales to smooth over (in pixels)
    if nargin < 2
        selectedLabels = [] ;
    end
    if nargin < 3
        sigmas = [5, 10, 20, 30, 50] ; 
    elseif isempty(sigmas)
        sigmas = [5, 10, 20, 30, 50] ; 
    end
    if nargin < 4
        steps = [1] ;
    elseif isempty(steps)
        steps = [1] ;
    end
    if nargin < 5 
        cdf_min = 0.01 ;    
        cdf_max = 0.999 ;
    elseif isempty(cdf_minmax)
        cdf_min = 0.01 ;    
        cdf_max = 0.999 ;
    else
        cdf_min = cdf_minmax(1) ;
        cdf_max = cdf_minmax(2) ;
    end
    if nargin < 6
        overwrite = false ;
    end
    
    if da.lookup.Count == 0 
        error('Before making gradients, define lookup map via dynamicAtlas.buildLookup()')
    end
    
    %% For each genoDir, build gradients 
    genoKeys = da.lookup.keys ;
    for geno_kk = 1:length(genoKeys)
        % Get dynamic atlas lookupMap for this genoDir
        daLUM = da.lookup(genoKeys{geno_kk}) ;
        
        %% For each labelDir, look at each embryo
        for key = daLUM.map.keys
            label = key{1} ;
            % Make an output directory for this label only in the genoDir
            outFigDir = fullfile(daLUM.genoDir, ...
                ['figures' filesep label filesep]) ;

            if ~exist(outFigDir, 'dir')
                mkdir(outFigDir)
            end

            % Select selected Channels: debug
            do_this_label = false ;
            if isempty(selectedLabels)
                do_this_label = true ;
            elseif any(contains(selectedLabels, label)) 
                do_this_label = true ;
            end
            
            if do_this_label
                lut = daLUM.map(label) ;
                disp('done building fileNames for this label')

                fileNames = lut.names ;
                fileDirs = lut.folders ;

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
                    % data = struct('im',[],'lim',[],'scale',[],'smooth',[]);

                    %% Load all files
                    % fbar = waitbar(0, '') ;
                    for k = 1:nFiles
                        disp(['Reading file ' num2str(k) '/' num2str(nFiles)]) 
                        % waitbar(k/nFiles, fbar, msg)
                        disp(fullfile(fileDirs{k}, fileNames{k}))
                        im = loadtiff(fullfile(fileDirs{k}, fileNames{k}));
                        disp('done loading')
                        
                        if size(im, 3) == 1
                            im = double(im) ;
                            im = im / max(im(:)) ;
                            disp('taking ecdf')
                            [f,x] = ecdf(im(:));
                            f1 = find(f>cdf_min, 1, 'first');
                            f2 = find(f<cdf_max, 1, 'last');
                            %imshow(tmp,[x(f1) x(f2)])
                            lim = [x(f1) x(f2)];
                            scale = mat2gray(im, double(lim));
                        else
                            disp('Scale each page of the image by its cdf')
                            im = double(im) / double(max(im(:))) ;
                            scale = im ;
                            for qq = 1:size(im, 3)
                                if mod(qq, 10) == 0                                   
                                   disp(['page' num2str(qq) ': taking cdf'])
                                end
                                page = double(squeeze(im(:, :, qq))) ;
                                page = page / max(page(:)) ;
                                [f,x] = ecdf(page(:));
                                f1 = find(f>cdf_min, 1, 'first');
                                f2 = find(f<cdf_max, 1, 'last');

                                %imshow(tmp,[x(f1) x(f2)])
                                lim = [x(f1) x(f2)];
                                scale(:, :, qq) = mat2gray(page, double(lim));
                            end
                        end
                        
                        %figure,imshow(scale)
                        disp('smoothing image')
                        smoothim = imfilter(scale, ...
                            fspecial('gaussian', kernel, sigma),...
                            'replicate');

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
                            if ~imsexist || overwrite
                                % Take gradients of the average of all data images
                                [X,Y] = meshgrid(1:size(smoothim, 1), 1:size(smoothim, 2));
                                %surf(X,Y,double(imavg)'), shading interp
                                [dX,dY] = gradient(double(permute(smoothim, [2, 1, 3])), step);
                                dX = dX / max(dX(:)) ;
                                dY = dY / max(dY(:));

                                gradmag = sqrt(dX.^2+dY.^2);
                                if view
                                    close all
                                    figure,
                                    %quiver(X(1:step:end,1:step:end),Y(1:step:end,1:step:end),dX,dY)
                                    pcolor(Y,X,gradmag), shading interp, axis equal, %colormap gray
                                    title([label ' gradient, $|\nabla I|$ '], 'Interpreter', 'Latex')
                                    pause(1)

                                    figure, 
                                    pcolor(Y,X,abs(dX)), shading interp, axis equal
                                    title([label ' $x$ gradient, $|\partial_x I|$ '], 'Interpreter', 'Latex')
                                    pause(1)

                                    figure, 
                                    pcolor(Y,X,abs(dY)), shading interp, axis equal
                                    title([label ' $y$ gradient, $|\partial_y I|$ '], 'Interpreter', 'Latex')
                                    pause(1)
                                end

                                % Save images of the gradients
                                disp('Writing smoothed im and gradients')
                                writeTiffPages(smoothim, sfn)
                                writeTiffPages(permute(gradmag, [2, 1, 3]), mfn) 
                                writeTiffPages(permute(abs(dX), [2, 1, 3]), xfn)
                                writeTiffPages(permute(abs(dY), [2, 1, 3]), yfn)
                                disp(['saved images for ' label ' ' num2str(k) '/' num2str(nFiles)])
                            end  
                        end
                        disp('Done with all steps')
                    end
                end
                disp(['Done with this sigma: ' sigma])
            end
        end
    end
