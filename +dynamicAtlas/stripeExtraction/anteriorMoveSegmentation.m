function [bw, bw2, placement] = ...
    anteriorMoveSegmentation(bw, placement, minsz, maxsz)
% Select a region of the BW segmented image that is anterior or posterior
% to current selection (indexed by 'placement') to find stripe 7
% segmentation.
%
% Parameters
% ----------
% bw : 
% placement : int
% minsz : float or int
% maxsz : float or int
%
% Returns
% -------
% button : waitforbuttonpress() output
%   what the current 
%
% NPMitchell 2020

% Enter loop where we are moving the chosen centroid
% anteriorly
while (strcmp(get(gcf, 'CurrentKey'), 'a') || ...
        strcmp(get(gcf, 'CurrentKey'), 'p'))
    
    % Move anteriorly or posteriorly
    if strcmp(get(gcf, 'CurrentKey'), 'a')
        placement = placement + 1 ;
    else
        placement = placement - 1 ;
    end
    
    % Move to anterior centroid by 1 object
    close all
    % Filter out small regions
    bw = bwareafilt(bw, [minsz maxsz]) ;
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
    bw2 = false(size(dcrop)) ;
    ind = min(length(sortind), ...
        max(1, length(sortind)-placement)) ;
    bw2(cc.PixelIdxList{sortind(ind)}) = true ;
    imshow(bw2)
    title(msg)
    disp(msg)
    disp('waiting for button press')
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
end