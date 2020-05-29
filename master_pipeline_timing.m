%% MASTER TIMING PIPELINE
% Master pipeline for building master timelines & timestamping fixed 
% samples using dynamic timeline.
% Demos of functionality are shown along the way.
%
% NPMitchell 2020

%% Add paths
% Add time_align_embryos directory to path so that dynamicAtlas package is
% available to use.
tlaDir = './';
addpath(tlaDir) ;
% Atlas path is where the dynamicAtlas resides
atlasPath = '/Users/npmitchell/Desktop/tmp/' ;

%% Build the dynamicAtlas
% Build dynamic atlas with all genotypes in the atlasPath
da = dynamicAtlas.dynamicAtlas(atlasPath) ;
% Or choose which genotypes to include in atlas (default=all of them)
genotypes = {'WT'} ;         
da = dynamicAtlas.dynamicAtlas(atlasPath, genotypes) ;

%% Build the lookupMaps for all genotypes
da.buildLookup() ;

%% Each map is now stored with a container
% example is
mapWT = da.lookup('WT') ;
% whose methods are
methods(mapWT) 
% for ex, to find embryos with a Runt stain
runts = mapWT.findLabel('Runt') 
% or equivalently
runts = mapWT.map('Runt') 
% To find embryos with a t=10 +/- 2 min
snaps = mapWT.findTime(10, 2) 
% To find Runt stains with a t=10 +/- 2 min
runtsnaps = mapWT.findLabelTime('Runt', 10, 2)

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

%% Create the gradients of pair rule genes for tensor analysis
% da.makeGradientImages()

% Note: Optionally, declare what scales to examine (in pixels)
% sigmas = [5, 10, 20, 30, 50] ;
% steps = [1] ;
% da.makeGradientImages(sigmas, steps)

%% Use dynamic datasets within the lookupMap to build master timeline
% To control how this is performed, toggle da.timeLineMethod
% Align dynamic runt nanobody data against each other 
da.makeMasterTimeline('WT', 'Runt')

%% Timestamp other data against the master timeline
% To control how this is performed, toggle da.timeStampMethod
% Here timestamp Runt fixed samples against runt timeline
% You can pass options through a struct if you wish.
Options = struct() ;
da.timeStamp('WT', 'Runt', Options)


