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
%   NPM added histequalize option on 27
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Name    = 'sqhCherryII_mean_rot_scaled_view1_intensity_scaled.tif';%'cylinder1_rot_scaled_view1.ome.tif';
Name = 'MAX_Cyl1_2_000000_c1_rot_scaled_view1.tif';


StackSize   = length(imfinfo(Name));
EdgeLength  = 15;   % Length of box edges in pixels; 
EdgeLength2 = 5;    % Length of boxedges in interpolated field
isf         = .4;   % image scaling factor. 
step        = 1;    % stp in timeframes.
smooth      = 1;    % set to 1 if gaussian smoothing is desired
KernelSize  = 10;    % Smoothing kernel size
sigma       = 2;  % standard deviation of gaussian kernel
histequilize= true ; % equilize histograms of raw data across each image (40x40 bins of image equilized)
im1 = imread(Name,1);
im1 = imresize(im1,isf,'bicubic');

%% 
% define the grid on which to compute the flow field
[X1,Y1] = meshgrid(EdgeLength/2:EdgeLength:size(im1,1)-EdgeLength/2,EdgeLength/2:EdgeLength:size(im1,2)-EdgeLength/2); 

for t = 1:  StackSize-step
    disp(['Running PIV on timestamp t = ' num2str(t)])

    % read the image and scale
    im1     = imread(Name,t);
    im2     = imread(Name,t+step);
    im1     = imresize(im1,isf,'bicubic'); % rescale image if desired
    im2     = imresize(im2,isf,'bicubic');
    %  im1 = im1(20:end-20,20:end-20);
    % im2 = im2(20:end-20,20:end-20);
   
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
    VX = VX-nanmean(VX(:));
    VY = VY-nanmean(VY(:));
    %VX = VX.*I;
    %VY = VY.*I;
    
    % Display image and overlay flow field.
%     imshow(im1',[])
%     hold on 
%     f = 10;
%     quiver(X1,Y1,f*VX,f*VY,0,'g-')
%     pause(.1)
    % records a movie
    %M(t)    = getframe(gcf); 
    if (~exist('PIV', 'dir'))
        mkdir('PIV')
    end
    save(sprintf('PIV/VeloT_%06d.mat',t),'VX','VY');

end
%play a movie;
%implay(M) 