function curv = extractStripeEdges(dat, maskfn, thres, minmaxsz, embryoID)

% Here, we cut the image in half, assuming that the time we are
% interested in finishes before stripe 7 extends halfway across
% the dorsal side
% Find the center pixel (rounded to nearest integer)
midx = round(0.5 * size(dat, 2)) ;
% Crop the image to right and center
dcrop = squeeze(dat(1, midx:end, :)) ;
dsz = [size(dat, 2), size(dat, 3)] ;
addx = midx ;
addy = 0 ;

% Unpack minsz, maxsz
minsz = minmaxsz(1) ;
maxsz = minmaxsz(2) ;

if exist(maskfn, 'file')
    bw2 = imread(maskfn);
else
    disp(['mask does not exist: ' maskfn])

    [bw2, bw] = autoStripe7(dcrop, thres, minsz, maxsz)  ;

    % Check if this automatic way is good enough
    move_on = false ;
    good_enough = false ;
    while ~move_on        
        % Restart the figure options
        close all
        imshow(bw2)
        msg = [embryoID ': Enter=OK, ', ...
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
                try
                    bw2 = BWmask .* bw ;
                catch
                    error('Previous mask is unusable')
                end
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
