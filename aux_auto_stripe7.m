% aux_auto_stripe7
% default stripe7 finding protocol
% 
% NPMitchell 2020
 
% % Find the center pixel (rounded)
% midx = round(0.5 * size(dat, 2)) ;
% % Crop the image to right and center
% dcrop = squeeze(dat(1, midx:end, :, tt)) ;
% dsz = [size(dat, 2), size(dat, 3)] ;
% addx = midx ;
% addy = 0 ;

% Obtain the back of stripe 7
% % Create BW image by thresholding
% datbw = false(size(dcrop)) ;
% datbw(dcrop > thres) = true ;
% posterior = bwareafilt(datbw, 1) ;
% bw = false(size(dcrop)) ;
% bw(posterior) = true ;
% % dilate and erode the binary image acting on the posterior reg
% se = strel('disk', 15) ;
% bw = imdilate(bw, se) ;
% bw2 = imerode(bw, se) ;
% % Now switch sign so posterior is false
% nbw = true(size(bw2)) ;
% nbw(bw) = false ;
% bw2 = nbw ;

% Obtain curve in front of stripe 7
datbw = false(size(dcrop)) ;
datbw(dcrop > thres) = true ;
% dilate and erode the binary image acting on the posterior reg
se2 = strel('disk', 2) ;
% se4 = strel('disk', 4) ;
bw = imdilate(datbw, se2) ;
bw = imerode(bw, se2) ;
% Filter out small regions
regions_exist = false ;
minsz_tmp = minsz ;
while ~regions_exist
    bwtest = bwareafilt(bw, [minsz_tmp maxsz]) ;
    if any(bwtest)
        regions_exist = true ;
    else
        minsz_tmp = minsz_tmp * 0.9 ;
    end
end
bw = bwtest ;

% Prepare for case where we use previous mask
prevmaskfn = fullfile(stripeDir, ['mask_' ...
             sprintf('%04d', tt-1) '.tif']); 

% Get mean of back object as first pass
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
placement = 0 ;  % seek maximum initiallye
bw2(cc.PixelIdxList{sortind(end-placement)}) = true ;
