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
tlaDir = '/Volumes/minimalData/code/';
% tlaDir = '/Users/mattlefebvre/Desktop/Code/code';
cd(fullfile(tlaDir)) ;
addpath(genpath('dynamicAtlasCode')) ;
cd('dynamicAtlasCode')
addpath(genpath('+dynamicAtlas'))
%import dynamicAtlas
% Atlas path is where the dynamicAtlas resides
% atlasPath = '/Volumes/minimalData/Atlas_Data' ;
atlasPath = '/Users/npmitchell/Desktop/Atlas_Data/' ;

%% Build the dynamicAtlas
% Build dynamic atlas with all genotypes in the atlasPath
% da = dynamicAtlas.dynamicAtlas(atlasPath) ;
% Or choose which genotypes to include in atlas (default=all of them)
da = dynamicAtlas.dynamicAtlas(atlasPath, {'WT'}) ;

%% Use dynamic datasets within the lookupMap to build master timeline
% To control how this is performed, toggle da.timeLineMethod
% Align dynamic runt nanobody data against each other 
Options.apCijFrac = 0.0 ;
Options = struct();
da.makeMasterTimeline('WT', 'Runt', Options)
%da.makeMasterTimeline('Tlrm9','Eve',Options)
