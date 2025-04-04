function makeGradientImages(da, selectedLabels, sigmas, steps, fixedOnly, ...
                            cdf_minmax, overwrite)
% MAKEGRADIENTIMAGES(da, selectedLabels, sigmas, steps,...
%                    cdf_minmax, overwrite)
%   Compute gradients of image intensities averaged over different 
% scales and store them on disk
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
% fixedOnly : boolean 
%   whether to make the gradient images on all (0) or fixed images only (1)
% steps : int array (optional, default=[1])
% cdf_minmax : length 2 float array (optional, default=[0.01, 0.999])
%   intensity cumulative distribution function limits for intensity
%   clipping each image
% overwrite : bool (otpional, default=false)
%   overwrite currently existing TIFF images
%
% Returns
% -------
% <nothing returned>
%
% Outputs
% -------
% 
%
% NPMitchell 2020
% Edited by Vishank Jain-Sharma 2025

% Declare what scales to smooth over (in pixels)
if nargin < 2
    selectedLabels = [] ;
end
if nargin < 3
    sigmas = [20] ; 
elseif isempty(sigmas)
    sigmas = [20] ; 
end
if nargin < 4
    steps = [1] ;
elseif isempty(steps)
    steps = [1] ;
end

%MODIFIED 2025/01/22 to add this variable in
if nargin < 5
    fixedOnly = 0 ;
elseif isempty(fixedOnly)
    fixedOnly = 0 ;
end

if nargin < 6 
    cdf_min = 0.01 ;    
    cdf_max = 0.999 ;
elseif isempty(cdf_minmax)
    cdf_min = 0.01 ;    
    cdf_max = 0.999 ;
else
    cdf_min = cdf_minmax(1) ;
    cdf_max = cdf_minmax(2) ;
end
if nargin < 7
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

            %% Now load all filenames and make gradients

            % Count the file names to do
            % for iii = 1:length(fileNames)
            %     disp(fileNames(iii))
            % end
            nFiles = length(fileNames);

            % Perform smoothing/gradients on each file
            for k = 1:nFiles

                %MODIFIED 2025/01/22: skip making gradients for live
                %datasets, if instructed to do so by fixedOnly argument
                if lut.nTimePoints(k) > 1
                    if fixedOnly
                        continue;
                    end
                end

                % Check if images are already done
                must_redo = false ;
                disp(['Checking if ' label ' ' lut.embryoIDs{k} ' gradients exist on disk'])   
                for ii = 1:length(sigmas)
                    sigma = sigmas(ii) ;
                    for jj = 1:length(steps)
                        step = steps(jj) ;
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
                            must_redo = true ;
                        end
                    end
                end

                % If we need to compute any images for the first
                % time, do so here
                if must_redo

                    % We must create new gradients/smoothed images, so
                    % load the TIFF file now
                    disp(['Reading file ', num2str(k), '/', ...
                            num2str(nFiles), ': ', lut.embryoIDs{k}]) 
                    disp(fullfile(fileDirs{k}, fileNames{k}))
                    im = loadtiff(fullfile(fileDirs{k}, fileNames{k}));
                    disp('done loading')

                    % Set intensity limits before smoothing/gradients
                    % 'scale' is the variable of the clipped image
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

                            %MODIFIED 2025/01/22 
                            %
                            %downsample the images before doing these
                            %operations, to save memory
                            ds_fac = 8;
                            scale_ds = 1 / ds_fac;
                            sz_1 = size(page,1);
                            sz_2 = size(page,2);
                            page_ds = imresize(page, [round(sz_1*scale_ds), round(sz_2*scale_ds)]);

                            %compute on downsampled version
                            [f,x] = ecdf(page_ds(:));
                            %[f,x] = ecdf(page(:)); %ORIGINAL

                            f1 = find(f>cdf_min, 1, 'first');
                            f2 = find(f<cdf_max, 1, 'last');

                            %imshow(tmp,[x(f1) x(f2)])
                            lim = [x(f1) x(f2)];

                            %compute on downsampled version
                            img_ds = mat2gray(page_ds, double(lim));
                            %resize to original size
                            img = imresize(img_ds,[sz_1,sz_2]);
                            scale(:, :, qq) = img;
                            %
                            %scale(:, :, qq) = mat2gray(page, double(lim)); %ORIGINAL
                            %
                            %%%%%
                        end
                    end

                    % Examine on all scales in kernel and gradient
                    for ii = 1:length(sigmas)
                        % Smooth the data with a gaussian filter
                        disp(['Smoothing with sigma = ' num2str(sigma)])
                        sigma = sigmas(ii) ;
                        kernel = 5*sigma;

                        % SECOND CHECK FOR FILE EXISTENCE 
                        % Before running imfiler, ask 'Does this sigma 
                        % file already exist for al smoothing/gradient 
                        % steps?'
                        must_redo_sigma = false ;
                        for jj = 1:length(steps)
                            step = steps(jj) ;
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
                                must_redo_sigma = true ;
                            end
                        end

                        if must_redo_sigma
                            %figure,imshow(scale)
                            disp('smoothing image')
                            smoothim = imfilter(scale, ...
                                fspecial('gaussian', kernel, sigma),...
                                'replicate');

                            %% Iterate over all gradient steps
                            for jj = 1:length(steps)
                                step = steps(jj) ;
                                disp([lut.embryoIDs{k} ' > Gradient with step ' num2str(step)])
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


                                    %%%%
                                    %
                                    %MODIFIED 2025/01/22
                                    %
                                    %Compute the gradient on a downsampled
                                    %version, then upsample
                                    ds_fac = 8;
                                    scale_ds = 1/ds_fac;
                                    step_ds = step * scale_ds;

                                    %size of the original image
                                    sz_1 = size(smoothim,1);
                                    sz_2 = size(smoothim,2);
                                    sz_3 = size(smoothim,3);
                                    
                                    %resize the image in the horizontal
                                    %dimensions
                                    
                                    %accounts for 3D (live) vs 2D (fixed)
                                    %stacks
                                    if (sz_3 > 1)
                                        smoothim_ds = imresize3(smoothim, [round(sz_1*scale_ds) round(sz_2*scale_ds) sz_3]);
                                    else
                                        smoothim_ds = imresize(smoothim, [round(sz_1*scale_ds) round(sz_2*scale_ds)]);
                                    end

                                    %defines permuted image explicitly
                                    smoothim_ds_permuted = double(permute(smoothim_ds, [2, 1, 3]));

                                    %defines gradient arrays explicitly
                                    %same size as the 
                                    dX_ds = 0 * double(smoothim_ds_permuted);
                                    dY_ds = 0 * double(smoothim_ds_permuted);

                                    %iterates over individual images to
                                    %compute the gradient on 2D   
                                    %
                                    %does so on downsampled images
                                    num_imgs = size(smoothim_ds_permuted,3);

                                    %successively computes the gradients
                                    %on the 2D slices, since only X and Y
                                    %components are important
                                    for n = 1 : num_imgs
                                        slice_n = squeeze(smoothim_ds_permuted(:,:,n));
                                        [dX_ds_n,dY_ds_n] = gradient(slice_n, step_ds);
                                        dX_ds(:,:,n) = dX_ds_n;
                                        dY_ds(:,:,n) = dY_ds_n;
                                    end

                                    %compute the gradient
                                    %[dX_ds,dY_ds] = gradient(double(permute(smoothim_ds, [2, 1, 3])), step_ds);
                                    
                                    if (sz_3 > 1)
                                        %account for the permutation above
                                        dX = imresize3(dX_ds, [sz_2 sz_1 sz_3]);
                                        dY = imresize3(dY_ds, [sz_2 sz_1 sz_3]);
                                    else
                                        dX = imresize(dX_ds, [sz_2 sz_1]);
                                        dY = imresize(dY_ds, [sz_2 sz_1]);
                                    end
                                    
                                    %ORIGINAL
                                    %[dX,dY] = gradient(double(permute(smoothim, [2, 1, 3])), step);
                                    %      
                                    %%%%

                                    gradmag = sqrt(dX.^2+dY.^2);
                                    denominator = max(gradmag(:)) ;
                                    gradmag = gradmag / denominator ;
                                    dX = dX / denominator ;
                                    dY = dY / denominator ;

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
                            disp(['Done with this sigma: ' sigma])
                        end
                    end
                else
                    disp(' -> already on disk')
                end 
            end
        end
    end
end
