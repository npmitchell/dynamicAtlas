function alignVelocityMagnitude(qs, speedCurv)


% First load all PIV results into qs
qs.getPIV() ;

% for qq = 1:length(qs.meta.names)
%     ddir = fullfile(qs.meta.folders{qq}, qs.meta.names{qq}) ;
%     for tidx = 1:(qs.meta.nTimePoints-1)
%         
%     end
% end

for qq = 1:length(qs.meta.names)
    speed = zeros(qs.meta.nTimePoints(qq) - 1, 1) ;
    for tidx = 1:(qs.meta.nTimePoints(qq) - 1)
        speed(tidx) = sqrt(mean(mean(qs.piv.vx{qq}(:, :, tidx).^2 + qs.piv.vy{qq}(:, :, tidx).^2))) ;
    end
    
    collapseCurves
end