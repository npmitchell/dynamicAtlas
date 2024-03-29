%% MASTER TIMING PIPELINE
% Template pipeline for computing lagrangian-advected patterns.
% Here we demonstrate by advecting Runt pattern along average PIV field
% built solely from sqh-mCherry and comparing to a later timepoint.
%
% NPMitchell 2021

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
% tlaDir = '/Volumes/minimalData/code/';
% tlaDir = '/Users/mattlefebvre/Desktop/Code/code';
tlaDir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/dynamicAtlas/code';
cd(fullfile(tlaDir)) ;
addpath(genpath('dynamicAtlas')) ;
cd('dynamicAtlas')
addpath(genpath('+dynamicAtlas'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define atlasPath to be where the dynamicAtlas resides (the parent
% dynamicAtlas directory, not the project directory '+dynamicAtlas')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%atlasPath = '/Volumes/minimalData/Atlas_Data' ;
%atlasPath = '/Users/mattlefebvre/Desktop/WT_data_server/'
atlasPath = '/Volumes/Elements/Atlas_Data' ;

%% Build the dynamicAtlas
% Build dynamic atlas with all genotypes in the atlasPath
% da = dynamicAtlas.dynamicAtlas(atlasPath) ;
% Or choose which genotypes to include in atlas (default=all of them)
options = struct() ;
options.labels = {'sqh-mCherry', 'Runt'} ;
da = dynamicAtlas.dynamicAtlas(atlasPath, {'WT'}, options) ;

%% Search for dynamic data 
genotype = 'WT' ;
% label = 'sqh-mCherry' ;
label = 'Runt' ;
qs = da.findDynamicGenotypeLabel(genotype, label) ;

qs = findEmbryo(da, embryoID) ;

%% Build average flow field at each point in time
pivStack = qs.buildPIVStack(da, genotype, label) ;

%% Create an ensemble reference Runt image and advect along pathlines
% Note there is a 12 minute difference between PIV timeline and Runt
% To find Runt stains with a t=10 +/- 2 min 
runtsnaps = da.findGenotypeLabelTime('WT', 'Runt', 12, 2) ;
im0 = runtsnaps.getMeanData() ;

Options = struct() ;
Options.outDir = outDir ;
% first timepoint of the desired image sequence
Options.minT = round(qs.getMinTime) ;
% last timepoint of the desired image sequence
Options.maxT = round(qs.getMaxTime) ;
% Whether to plot a video of some Lagrangian tracer beads first
Options.plot_scatterpaths = false ;
% Apply optical flow to reference image im0
applyOpticalFlow(pivStack, im0, Options)

% Compare to mean over time
dmyk = 1 ;
for tt = round(qs.getMinTime):round(qs.getMaxTime)
    imfn = fullfile(outDir, 'meanRunt', [sprintf('%03d', dmyk), '.png']) ;
    if ~exist(imfn, 'file')
        runtsnaps = da.findGenotypeLabelTime('WT', 'Runt', 12+tt, 1) ;
        im0 = runtsnaps.getMeanData() ;
        if ~isempty(im0)
            imwrite(im0, imfn)
        end
    end
    dmyk = dmyk + 1 ;
end


%% OTHER FUNCTIONALITY -- Apply average flow field to a single Runt image
% Create a reference image to advect 
imfn = '/Volumes/Elements/Atlas_Data/WT/Runt/202001150004/MAX_Cyl1_2_000000_c1_rot_scaled_view1_ss04.tif' ;
t0_furrow = 12 ;
im0 = imread(imfn, t0_furrow);

outDir = fullfile('/Users/npmitchell/Desktop/tmp/WT_pullbackPathlines') ;
if ~exist(outDir, 'dir')
    mkdir(outDir)
end

Options = struct() ;
Options.outDir = outDir ;
% first timepoint of the desired image sequence
Options.minT = round(qs.getMinTime) ;
% last timepoint of the desired image sequence
Options.maxT = round(qs.getMaxTime) ;
% Whether to plot a video of some Lagrangian tracer beads first
Options.plot_scatterpaths = true ;
% Apply optical flow to reference image im0
applyOpticalFlow(pivStack, im0, Options)


%% Search for dynamic data with PIV
% qs = findFlowGenotypeLabel(da, genotype, label) ;
% da.makeMasterFlowField('WT', 'Runt', Options)


