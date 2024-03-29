function [bw, bw2, good_enough, move_on, placement] =...
    morphologicalOpLoop(bw, dcrop_orig, placement, thres, minmaxsz, msg)
%MORPHOLOGICALOPLOOP(bw, dcrop_orig, placement, thres, minmaxsz)
%   Loop to extract mask segmentation by iterative dilation & erosion
% 
% Parameters
% ----------
% bw : NxM binary array
%   pre-size-filtered binary mask
% dcrop_orig : original intensity data
%   original intensity data before thresholding to create mask
% placement : int
%   the rank AP position of the segmented region (stripe7, for ex) that is
%   true in the mask     
% thres : numeric threshold
%   threshold for binarizing the dcrop_orig data 
% minmaxsz : [minsz maxsz]
%   the min and max allowed size of a true region in the mask, in pixels
% 
% Returns
% -------
% bw : NxM binary array
%   pre-size-filtered binary mask
% bw2 : NxM binary array
%   size-filtered binary mask
% good_enough : bool
%   segmentation is good enough to save, should be true when returned
% move_on : bool
%   whether we are ready to continue with the pipeline
% placement : int
%   the rank AP position of the segmented region (stripe7, for ex) that is
%   true in the mask 
%
% NPMitchell 2020

% Enter loop where we continually morphologically
% refine by eroding + dilating
se2 = strel('disk', 2) ;
move_on = false ;
% Initially remember current state as input state
bw2 = bw ;
while (strcmp(get(gcf, 'CurrentKey'), 'e') || ...
        strcmp(get(gcf, 'CurrentKey'), 'd'))
    % Morphological operation
    if strcmp(get(gcf, 'CurrentKey'), 'e')
        bw = imerode(bw, se2) ;
    elseif strcmp(get(gcf, 'CurrentKey'), 'd')
        bw = imdilate(bw, se2) ;
    end
    close all

    % recall current state later as bwprev
    bwprevA = bw2 ;
    % Filter out small regions
    bw = bwareafilt(bw, minmaxsz) ;
    cc = bwconncomp(bw) ;
    rp = regionprops(cc) ;
    centry = zeros(length(rp), 1) ;
    for qq = 1:length(rp)
        centry(qq) = rp(qq).Centroid(2) ;
        % hold on;
        % plot(rp(qq).Centroid(1), rp(qq).Centroid(2), 'o')
    end
    % [~, ind] = max(centry) ;
    [~, sortind] = sort(centry) ;
    if placement > (length(sortind) -1)
        placement = length(sortind) - 1 ;
    elseif placement < 0
        placement = 0 ;
    end
    if isempty(sortind)
        disp('No regions found, resetting bw2=bw')
        datbw = false(size(dcrop_orig)) ;
        datbw(dcrop > thres) = true ;
        % dilate and erode the binary image acting on the posterior reg
        se2 = strel('disk', 2) ;
        % se4 = strel('disk', 4) ;
        bw = imdilate(datbw, se2) ;
        bw = imerode(bw, se2) ;
        % Filter out small regions
        bw2 = bwareafilt(bw, minmaxsz) ;
    else
        bw2 = false(size(bw)) ;
        bw2(cc.PixelIdxList{sortind(end-placement)}) = true ;
    end
    % Show the image
    subplot(2, 1, 1)
    imshow(bw2)
    subplot(2, 1, 2)
    imagesc(bw2 - bwprevA)
    colorbar()
    title(msg)
    disp(msg)
    disp('waiting for button press')
    button = waitforbuttonpress() ;
    good_enough = false ;
    if button && strcmp(get(gcf, 'CurrentKey'), 'return')
        disp('Good enough using masking of other centroid objects, saving mask')
        disp(['Saving mask to ' maskfn])
        imwrite(bw2, maskfn) ;
        move_on = true ;
    elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
        disp('Define manual ROI poly')
        move_on = true ;
    end
end