function segs = cutCurveIntoLeadingTrailingSegments(stripe)
% cutCurveIntoLeadingTrailingSegments()
% cut curve into left and right segments (leading and trailing segments)
%
% Assume clockwise motion of the curve
%
% NPMitchell 2020

[maxval, maxi] = max(stripe(:, 2)) ;
aa = stripe(1:maxi(1)+1, :) ;
bb = stripe(maxi(end):end, :) ;

% So now we have two parts, aa and bb
% Grab starting point for aa

[minval, mini] = min(stripe(:, 2)) ;
if mini(1) < maxi(1)
    % This is the easy case. 
    % Postpend the start of the curve to the end of bb
    bb = cat(1, bb, aa(1:mini(1), :)) ;
    % Now keep only trailing part of bb
    mxi = find(bb(:, 2) == maxval) ;
    mni = find(bb(:, 2) == minval) ;
    bb = bb(mxi(end):mni(1)-1, :) ;
elseif mini == maxi
    error('here')
else
    error('min == max for this curve. Handle here.')
end

segs{1} = aa ;
segs{2} = bb ;

% check it
% plot(aa(:,1), aa(:,2), '.') ;
% hold on;
% plot(bb(:,1), bb(:,2), '.') ;
