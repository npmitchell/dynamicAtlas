%% MASTER TIMING PIPELINE
% Master pipeline for building master timelines & timestamping fixed 
% samples using dynamic timeline.
% Demos of functionality are shown along the way.
%
% NPMitchell 2020

%% First mount the server onto your machine 
% for example, on a Mac: Apple+K afp://flydrive.synology.me
% mount minimalData/ 

%% Let's clear our environment
clc
clear
close all

%% Add paths (this part can be slow)
% Add time_align_embryos directory to path so that dynamicAtlas package is
% available to use.
tlaDir = './';
cd(fullfile(tlaDir)) ;
addpath(genpath('./')) ;
addpath(genpath('+dynamicAtlas'))
%import dynamicAtlas
% Atlas path is where the dynamicAtlas resides
% atlasPath = '/Volumes/minimalData/Atlas_Data/' ;
atlasPath = '/Volumes/Elements/Atlas_Data/' ;
% atlasPath = '/Users/mattlefebvre/Desktop/WT_data_server/'
%% Build the dynamicAtlas
% Build dynamic atlas with all genotypes in the atlasPath
% da = dynamicAtlas.dynamicAtlas(atlasPath) ;
% Or choose which genotypes to include in atlas (default=all of them)
da = dynamicAtlas.dynamicAtlas(atlasPath, {'WT'}) ;

sqh = da.findGenotypeLabel('WT','sqh-mCherry')  

%% The properties of da are:
properties(da)

%% Look at the methods of da
methods(da)

%% We can change which genotypes are included at any time via
% da.buildLookup({'WT', 'Eve'}) 

%% We can grab all the data from a label within a genotype 
qs = da.findGenotypeLabel('WT', 'Runt')
% Now we can see what data lies in this slice via
qs.meta
% and we can extract that data via 
qs.getData()
% The data is stored in a property
qsdata = qs.data ;

%% We can grab all the data associated with a given embryo
embryoID = '201905091604' ;
qs2 = da.findEmbryo(embryoID) ;
% Now qs is a queriedSample object with methods
qs.getData() ;
qs.getMeanData() ;

%% Each genotype's map is now stored with a container called da.lookup
% example is
mapWT = da.lookup('WT') ;
% whose methods are
methods(mapWT) 
% To find embryos with a Runt stain, for ex,
runts = mapWT.findLabel('Runt')  
runts = da.findGenotypeLabel('WT','Runt')  
% To find embryos with a t=10 +/- 2 min
snaps = mapWT.findTime(30, 5) 
snaps = da.findGenotypeLabelTime('WT', 'Runt', 10, 2)
% To find Runt stains with a t=10 +/- 2 min
runtsnaps = mapWT.findLabelTime('Runt', 20, 2)


timeseq = [] ;
for tt = 1:4
    snaps = mapWT.findTime(tt, 0.5) ; 
    timeseq(tt, :, :) = snaps.getMeanData() ;
end

% Note: mapWT has a property called map that contains a key for each label
% For ex, 'Runt' is a key that returns a struct with fields
%   times : 1xN float array
%       the timestamps for each embryo
%   uncs : 1xN float array
%       the uncertainty in time for each embryo
%   folders : 1xN string cell
%       the folder containing the pullback for each embryo
%   names : 1xN string cell
%       file name within folder of the pullback image
%   embryoIDs : 1xN string cell
%       unique date identifier for each embryo
%   nTimePoints : 1xN int array
%       number of timepoints in the pullback

%% The lookup maps can be queried also.
% For ex, da.lookup('WT') is itself a map
da.lookup('WT')
da.lookup('WT').findTime(10, 2)


%% Get metric of the embryo mesh and apply to the PIV fields
rect = read_ply('rect.ply') ;
erec = read_ply('embryo_rect.ply') ;
rcut = read_ply('embryo_rect_noglue.ply') ;
ecoarse = read_ply('embryo_coarse.ply') ;
ecut = read_ply('coarse_cylindrical_chart_ventral_cut.ply') ;

[gg, bb] = constructFundamentalForms(rcut.f, rcut.v, rect.v(:, 1:2)) ;
% rect is a rectilinear planar triangulation of the embryo surface as
% topological disk (shared edge along DV direction)
% embryo_rect is the embedding of rect as a glued topological cylinder
% embryo_rect_noglue has a cut for the seam in DV
% embryo_coarse is a MeshLab reconstruction in embedding space (3D)
% coarse_cylindrical_chart_ventral_cut is a mapping of embryo_coarse to the
% plane with a cut on the ventral side (dorsal is on the midline)
save('embryo_rect_fundamentalForms.mat', 'gg', 'bb')

%% Now construct metric in PIV coord grid
% standard size for pullbacks in atlas: 696 x 820 

% Built-in pullback image size, hard-coded here:
% ----------------------------------------------
% Original image is 2000x2000, then resized is 0.4*(Lx,Ly), then PIV is 15x
% smaller than the resized image. 
% fullsizeCoordSys: 1738x2050
% resizeCoordSys: 696 x 820
% pivCoordSys: 46 x 54
EdgeLength = 15;
% rescaleFactor = 0.4;
% szX_orig = 1738 ;
% szY_orig = 2050 ;

% resized dimensions of piv grid --> nearly (szX_orig, szY_orig) * isf
rescaleFactor = 0.4 ;
szX = 696 ;
szY = 820 ;

% PIV evaluation coordinates in resized pixels
[X0,Y0] = meshgrid(EdgeLength/2:EdgeLength:(szX-EdgeLength/2), ...
    EdgeLength/2:EdgeLength:(szY-EdgeLength/2)); 

% Resize embryo rect to same size as downsampled images for PIV
rectPIVscale = rect ;
rectPIVscale.v(:, 1) = rect.v(:, 1) * 696/218 ;
rectPIVscale.v(:, 2) = rect.v(:, 2) * 820/257 ;
plywrite('./rect_PIVImageScale.ply', rectPIVscale.f, rectPIVscale.v)

[gg, bb] = constructFundamentalForms(rcut.f, rcut.v, rectPIVscale.v(:, 1:2)) ;

% What is the metric at the sites of PIV evaluation?
TRpiv = triangulation(rectPIVscale.f, rectPIVscale.v(:, 1:2)) ;
% find PIV evaluation coordinate faces in TRpiv
faceID = pointLocation(TRpiv, [X0(:), Y0(:)] ) ;
ggPIV = gg(faceID) ;
bbPIV = bb(faceID) ;
readme = ['ggPIV and bbPIV are cell arrays of the first and second ', ...
    'fundamental form evaluated at the PIV coordinate locations in a ', ...
    'resized image (0.4 * original pixel dimensions) when compared to ', ...
    'the embedding space of the embryo surface. The ith cell is the ', ...
    'fundamental form of [X0(i), Y0(i)] evaluation point -- note that X0 and Y0 are unravelled for this indexing. \n' ] ;
readme = [readme 'X0 and Y0 are the resized image pixel positions of the PIV evaluation points. \n'] ;
readme = [readme 'gg and bb are the first and second fundamental forms ', ...
    'for each face of rect_embryo_noglue.ply compared to the pullback rect_PIVImageScale.ply'] ;

save('embryo_rectPIVscale_fundamentalForms.mat', ...
    'gg', 'bb', 'ggPIV', 'bbPIV', 'X0', 'Y0', 'readme')

% save('./rect_




%%
trisurf(triangulation(ecut.f, ecut.v)) ;
view(2)
axis equal
