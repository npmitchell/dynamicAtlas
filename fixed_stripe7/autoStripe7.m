function [bw2, bw] = autoStripe7(dcrop, thres, minsz, maxsz) 
% autoStripe7()
% default stripe7 finding protocol
% 
% NPMitchell 2020
 
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
dmyk = 0 ;
while ~regions_exist
    disp(['checking that regions exist: pass ' num2str(dmyk)])
    bwtest = bwareafilt(bw, [minsz_tmp maxsz]) ;
    if any(bwtest(:))
        regions_exist = true ;
    else
        % get largest region size and display it
        regs = regionprops(bw) ;
        regsizes = zeros(length(regs), 1) ;
        for dmyq = 1:length(regs)
            regsizes(dmyq) = regs(dmyq).Area ;
        end
        biggestreg = max(regsizes) ;
        disp(['Largest region size is ' num2str(biggestreg)])
        disp('Using that size as minsize')
        minsz_tmp = biggestreg - 1 ;
    end
    dmyk = dmyk + 1 ;
end
bw = bwtest ;

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
