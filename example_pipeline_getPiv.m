%% Let's clear our environment
clc
clear
close all

%% Add paths (this part can be slow)
% Add time_align_embryos directory to path so that dynamicAtlas package is
% available to use.
tlaDir = '/path/to/dynamicAtlas/githubRepo/';
% tlaDir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/dynamicAtlas/code';

cd(fullfile(tlaDir)) ;
addpath(genpath('dynamicAtlas')) ;
cd('dynamicAtlas')
addpath(genpath('+dynamicAtlas'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define atlasPath to be where the dynamicAtlas resides (the parent
% dynamicAtlas directory, not the project directory '+dynamicAtlas')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% atlasPath = '/Volumes/minimalData/Atlas_Data' ;
atlasPath = '/Users/npmitchell/Desktop/Atlas_Data/';

%% Build the dynamicAtlas
% Build dynamic atlas with just WAT genotype in the atlasPath
options = struct() ;
options.labels = { 'tmp-CAAX-mCherry-on3'} ;  
options.prepend = 'MAX_Cyl1_2_00000*_c*_rot_scaled_view1_EVERY_FOUR' ; 
da = dynamicAtlas.dynamicAtlas(atlasPath, {'WT'}, options) ;

%% Build flow field at each point in time using PIVLab
qs = da.findEmbryo('201807021620') ;
options = struct() ;
options.method = 'pivlab'; % 'default' is the other option
qs.ensurePIV(options) ;