function energy = springEnergy1D(xx, BL, fixed_ind, fixed_x) 
% energy = springEnergy1D(xx, BL, fixed_ind, fixed_x) 
%
% Compute spring energy as distance between positons of nodes that share
% bonds, minus rest bond length, squared. Here we subtract the
% bondLength BEFORE squaring, so bonds are directed and have an oriented 
% bond length!
%
% NOTE: xx has only mobile dof 
% Insert fixed positions into xx
%
% Parameters
% ----------
% xx : (#nodes-Q) x 1 numeric, where Q=#immobile nodes
%   current positions of the nodes that are connected by bonds
% BL : #bonds x 2 int
%   bond connections between nodes, such that BL(ii, 1) is connected to
%   BL(ii, 2) with signed rest length bondLength0 (ie at rest, 
%   xx(BL(ii,2)) == xx(BL(ii, 1) + bondLength0(ii) ; 
%   The strength of bond ii is given by stiffness(ii).
% fixed_ind : Qx1 int
%   indices of nodes/vertices that are considered immobile
% fixed_x : Qx1 numeric
%   positions of the indices which are not dof
%
% Returns
% -------
% energy : 
%   1/2 * sum(dL.^2), where dL is the change in bond length from the
%   rest bond length bondLength0. Fixed nodes can contribute energy but
%   will not be moved in the minimization of this functional. For ex,
%     options = optimset('PlotFcns',@optimplotfval, 'TolX', 1e-4, 'TolFun', 1e-6);
%     x0 = i_tau0j(:, 2) ;
%     % pop the indices of the fixed times from the array x0
%     fixed_ind = find(i_tau0j(:, 1) == hard) ;
%     fixed_xx = i_tau0j(fixed_ind, 2) ;
%     x0(fixed_ind) = [] ;
%     fun = @(x)springEnergy1D(x, BL, fixed_ind, fixed_xx);
%     xf = fminsearch(fun, x0, options) ;
%
%
xnew1 = xx(1:fixed_ind(1)-1) ;
xnew2 = xx(fixed_ind(1):end) ;
xnew = [xnew1; fixed_x; xnew2 ] ;
energy = 0.5 * sum(abs(xnew(BL(:, 2)) - xnew(BL(:, 1))).^2) ;

