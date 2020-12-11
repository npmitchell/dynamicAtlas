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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Name    = 'sqhCherryII_mean_rot_scaled_view1_intensity_scaled.tif';%'cylinder1_rot_scaled_view1.ome.tif';
Name = 'CAAX-mCherry_singlelayer_rot_scaled_view1_intensity_scaled.tif';


StackSize   = length(imfinfo(Name));
EdgeLength  = 15;   % Length of box edges in pixels; 
EdgeLength2 = 5;    % Length of boxedges in interpolated field
isf         = .4;   % image scaling factor. 
step        = 1;    % step in timeframes.
smooth      = 1;    % set to 1 if gaussian smoothing is desired
KernelSize  = 10;    % Smoothing kernel size
sigma       = 2;  % standard deviation of gaussian kernel
im1 = imread(Name,1);
im1 = imresize(im1,isf,'bicubic');
%im1 = im1(20:end-20,20:end-20);
%% 
% define the grid on which to compute the flow field
[X1,Y1] = meshgrid(EdgeLength/2:EdgeLength:size(im1,1)-EdgeLength/2,EdgeLength/2:EdgeLength:size(im1,2)-EdgeLength/2); 

for t = 1:  StackSize-step

    % read the image and scale
    im1     = imread(Name,t);
    im2     = imread(Name,t+step);
    im1     = imresize(im1,isf,'bicubic'); % rescale image if desired
    im2     = imresize(im2,isf,'bicubic');
  %  im1 = im1(20:end-20,20:end-20);
   % im2 = im2(20:end-20,20:end-20);
    % compute the piv flow field
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
    imshow(im1',[])
    hold on 
    f = 10;
    quiver(X1,Y1,f*VX,f*VY,0,'g-')
    save(sprintf('PIV/VeloT_%06d.mat',t),'VX','VY');
    pause(.1)
    % records a movie
    M(t)    = getframe(gcf); 
end
%play a movie;
implay(M) 