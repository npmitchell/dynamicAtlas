function energy = springEnergy1D(xx, BL, fixed_ind, fixed_x) 
% xx has only mobile dof 
% Insert fixed positions into xx
xnew1 = xx(1:fixed_ind(1)-1) ;
xnew2 = xx(fixed_ind(1):end) ;
xnew = [xnew1; fixed_x; xnew2 ] ;
energy = 0.5 * sum(xnew(BL(:, 1)) - xnew(BL(:, 2))).^2 ;