function PIVTimeseries(inputDir, options) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Using GetPIV we computes the flow field estimate based on an image sequence 
%   im1 & im2 are assumed to have the same dimensions. The grid X1,Y1
%   is assmued to be contained in the image domain with finite
%   EdgeLength defining size of PIV box. 
%   This script generates a movie displaying the original image with flow
%   field overlayed. 
%   
%   Written by: Sebastian J Streichan, KITP, February 14, 2013
%   NPM added histequalize option and made function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

method = 'default' ; % 'pivlab' or 'default'
%Name    = 'sqhCherryII_mean_rot_scaled_view1_intensity_scaled.tif';%'cylinder1_rot_scaled_view1.ome.tif';
Name = 'MAX_Cyl1_2_000000_c1_rot_scaled_view1.tif';
EdgeLength  = 15;   % Length of box edges in pixels; 
% EdgeLength2 = 5;    % Length of boxedges in interpolated field
isf         = .4;   % image scaling factor. 
step        = 1;    % stp in timeframes.
smooth      = 1;    % set to 1 if gaussian smoothing is desired
KernelSize  = 10;    % Smoothing kernel size
sigma       = 2;  % standard deviation of gaussian kernel
histequilize= true ; % equilize histograms of raw data across each image (40x40 bins of image equilized)

if nargin > 1
    if isfield(options, 'method')
        method = options.method ;
    end
    if isfield(options, 'Name')
        Name = options.Name ;
    end
    if isfield(options, 'EdgeLength')
        EdgeLength = options.EdgeLength ;
    end
    if isfield(options, 'isf')
        isf = options.isf ;
    end
    if isfield(options, 'smooth')
        smooth = options.smooth ;
    end
    if isfield(options, 'KernelSize')
        KernelSize = options.KernelSize ;
    end
    if isfield(options, 'sigma')
        sigma = options.sigma ;
    end
    if isfield(options, 'histequilize')
        histequilize = options.histequilize ;
    end
end

% Now that name is defined, find stacksize
StackSize   = length(imfinfo(fullfile(inputDir, Name)));

%% Define the grid on which to compute the flow field
if strcmpi(method, 'default')
    % Get image size for defining grid on which to evaluate PIV
    im1 = imread(fullfile(inputDir, Name), 1);
    im1 = imresize(im1,isf,'bicubic');

    % define the grid on which to compute the flow field
    [X1,Y1] = meshgrid(EdgeLength/2:EdgeLength:size(im1,1)-EdgeLength/2,EdgeLength/2:EdgeLength:size(im1,2)-EdgeLength/2); 
end

%% Get full data for page-by-page analysis
disp(['Reading tiff: ' fullfile(inputDir, Name)])
dat = loadtiff(fullfile(inputDir, Name)) ;

%% Compute PIV for each frame
dt = dlmread(fullfile(inputDir, 'dt.txt')) ;
pix2um = dlmread(fullfile(inputDir, 'pix2um.txt')) ;
for t = 1:  StackSize-step
    disp(['Running PIV on timestamp t = ' num2str(t)])

    % read the image and scale
    im1     = dat(:,:,t);
    im2     = dat(:,:,t+step);
    
    if isf ~= 1.0
        im1     = imresize(im1,isf,'bicubic'); % rescale image if desired
        im2     = imresize(im2,isf,'bicubic');
        % im1 = im1(20:end-20,20:end-20);
        % im2 = im2(20:end-20,20:end-20);
    end
    
    if strcmpi(method, 'default')
       
       if histequilize
           im1 = histeq(im1) ;
           im2 = histeq(im2) ;
       end

        % compute the piv flow field
        % If GetPIV is not found, then run this line:
        % >> addpath addpath '/Volumes/minimalData/code/dynamicAtlasCode/+dynamicAtlas/piv_handling'
        [VX,VY] = GetPIV(im1,im2,X1,Y1,EdgeLength); 

        % smooth if desired
        if smooth == 1
            VX  = imfilter(VX,fspecial('gaussian',KernelSize,sigma))/step;
            VY  = imfilter(VY,fspecial('gaussian',KernelSize,sigma))/step;
        end
        % VX = VX-nanmean(VX(:));
        % VY = VY-nanmean(VY(:));
        % VX = VX.*I;
        % VY = VY.*I;

        % Inspect velocity
        clf
        subplot(1, 2, 1)
        imagesc(VX)
        colorbar; title(['$v_x$, t=' num2str(t)], 'interpreter', 'latex')
        subplot(1, 2, 2)
        imagesc(VY)
        colorbar; title(['$v_y$, t=' num2str(t)], 'interpreter', 'latex')
        pause(0.1)

        % Display image and overlay flow field.
        %     imshow(im1',[])
        %     hold on 
        %     f = 10;
        %     quiver(X1,Y1,f*VX,f*VY,0,'g-')
        %     pause(.1)
        % records a movie
        %M(t)    = getframe(gcf); 
        if ~exist(fullfile(inputDir, 'PIV'), 'dir')
            mkdir(fullfile(inputDir, 'PIV'))
        end
        save(sprintf(fullfile(inputDir, 'PIV', 'VeloT_%06d.mat'),t),'VX','VY');

    elseif strcmpi(method, 'pivlab')
        
        optfn = fullfile(inputDir, 'PIVlab_settings.mat') ;
        if exist(optfn, 'file')
           opts = load(optfn) ;
        else
            opts = struct() ;
            
            % Standard PIV Settings
            opts.intArea1        = EdgeLength*8 ;
            opts.step            = round(opts.intArea1 * 0.5) ;
            opts.subpixFindr     = 1  ;
            opts.mask            = [] ;
            opts.roi             = [] ;
            opts.numPasses       = 4  ;
            opts.intArea2        = EdgeLength*4 ;
            opts.intArea3        = EdgeLength*2 ;
            opts.intArea4        = EdgeLength ;
            opts.repeat          = 1  ;
            opts.disAuto         = 0  ;
            % Image proc
            opts.clahe = histequilize ;
            opts.claheW          = 40 ;
            opts.highPass        = 0  ;
            opts.highPassSz      = 15 ;
            opts.clipping        = 0  ;
            opts.wiener          = 0  ;
            opts.wienerW         = 3  ;
            % Post-processing        
            opts.calu            = 1. ;
            opts.calv            = 1  ;
            opts.valid_vel       = [] ;
            opts.do_stdev_check  = 1 ;
            opts.stdthresh       = 5  ;
            opts.do_local_median = 0  ;
            opts.neigh_thresh    = 2  ;
        end
        
        [VX, VY, VX_filt, VY_filt, X0, Y0] = ...
            getPIVLab(im1, im2, opts) ;
        
        % Inspect velocity
        clf  
        quiver(X0, Y0, 10*VX, 10*VY, 0)
        axis equal
        title(['$\vec{v}$, t=' num2str(t)], 'interpreter', 'latex')
        xlim([0, 696]) ;
        ylim([0, 820]) ;
        pause(0.0001)

        % Display image and overlay flow field.
        %  imshow(im1',[])
        %  hold on 
        %  f = 10;
        %  quiver(X0,Y0,f*VX,f*VY,0,'g-')
        %  pause(.1)
        % records a movie
        %M(t)    = getframe(gcf); 
        if ~exist(fullfile(inputDir, 'PIVlab'), 'dir')
            mkdir(fullfile(inputDir, 'PIVlab'))
        end
        fn = sprintf(fullfile(inputDir, 'PIVlab', 'VeloT_fine_%06d.mat'),t) ; 
        VX = VX' ;
        VY = VY' ;
        VX_filt = VX_filt' ;
        VY_filt = VY_filt' ;
        X0 = X0' ;
        Y0 = Y0' ;
        convert_to_um_per_min = pix2um / (isf * dt) ;
        save(fn,'VX','VY','VX_filt','VY_filt','X0','Y0', 'convert_to_um_per_min');
        
        % Save the settings used for piv as a settings.mat
        if ~exist(optfn, 'file')
            save(optfn, 'opts') ;
        end
    end
end
%play a movie;
%implay(M) 