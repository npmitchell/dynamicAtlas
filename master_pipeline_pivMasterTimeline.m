%% Script for alignment of dynamic datasets based on displacement


% tlaDir = '/Users/mattlefebvre/Desktop/Code/code';
tlaDir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/dynamicAtlas/code';
% tlaDir = '/mnt/data/code/';

cd(fullfile(tlaDir)) ;
addpath(genpath('dynamicAtlas')) ;
cd('dynamicAtlas')
addpath(genpath('+dynamicAtlas'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define atlasPath to be where the dynamicAtlas resides (the parent
% dynamicAtlas directory, not the project directory '+dynamicAtlas')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% atlasPath = '/Volumes/minimalData/Atlas_Data' ;
atlasPath = '/Volumes/DOGIC/Atlas_Data' ;
atlasPath = '~/Desktop/Atlas_Data' ;

options = struct() ;
genotype = 'twistPlus' ;
label = 'sqhGFP' ;
options.labels = { label } ; 
options.timeStampMethod = 'pullbackPathlineDisplacement' ;
da = dynamicAtlas.dynamicAtlas(atlasPath, {genotype}, options) ;
da.makeMasterTimeLinePIV('twistPlus', 'sqhGFP')

%% Get snapshots of piv for figures at three different times
qs = da.findDynamicGenotypeLabel('twistPlus', 'sqhGFP') ;
piv = qs.getPIV(struct('method', 'pivlab')) ;
exptIDs = qs.meta.embryoIDs ;


for tt = [15, 25, 35] 
    qsT = da.findDynamicGenotypeLabelTime('twistPlus', 'sqhGFP', tt, 1) ;
    options = struct('method', 'pivlab') ;
    meanPIV = qsT.getMeanPIV(options) ;

    clf
    vx = meanPIV.vx * 0.2619 ;
    vy = meanPIV.vy * 0.2619 ;
    X0 = meanPIV.X0 * 0.2619 ;
    Y0 = meanPIV.Y0 * 0.2619 ;
    speed = sqrt(vx.^2 + vy.^2) ;
    % subsample for quiver
    
    fac = 0.1 ;
    X0sub = imresize(X0, fac) ;
    Y0sub = imresize(Y0, fac) ;
    vxsub = imresize(vx, fac) ;
    vysub = imresize(vy, fac) ;    
    
    figure('Units', 'centimeters', 'position', [0,0,6,6])
    imagesc(unique(X0), unique(Y0), speed') ;
    hold on;
    axis equal
    quiver(X0sub, Y0sub, vxsub * 10, vysub * 10, 0)
    axis equal 
    axis tight
    cb = colorbar() ;
    ylabel(cb, 'speed [\mum/min]')
    axis off
    pause(1)
end
