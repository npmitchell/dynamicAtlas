%% TIMING SCRIPT -- from linear offset of live data wrt master timeline
% Adding timestamps to live data 
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
cd(fullfile(tlaDir)) ;
disp(['Adding path to dynamicAtlasCode'])
addpath(genpath('dynamicAtlasCode')) ;
cd('dynamicAtlasCode')
disp(['Adding path to +dynamicAtlasCode'])
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


%% Some offsets computed via |piv|(t) curve collapse
% Hannah's twist control shh::GFP movies:

sqhT = struct('202006261115', -1, ...
 '202007011145', 21, ...
 '202007081130', 13, ...
 '202007091200', 21, ...
 '202008131005', 10, ...
 '202009041145', 11) ;

% Runt::LlamaTag-GFP movies:

runtT = struct('202001210044', 6, ...
 '202001142033', 27, ...
 '202001150004', 18, ... % # Noah's master timeline embryo
 '202001141730', 13) ;

% Runt movies 202001141943, 202001210000 start too late for the inflection point-based PIV method to work. But since you have relative timing between them and the Runt master timeline based on stripes, this should be no problem.

% Eve::YFP movies:

eveT = struct('201406122111', 23, ...
  '202001201216', -1, ...
  '202001201716', 24, ...
  '202001202210', 10) ;


%% The properties of da are:
properties(da)

%% Look at the methods of da
methods(da)

%% We can change which genotypes are included at any time via
% da.buildLookup({'WT', 'Eve'}) 

%% We can grab all the data from a label within a genotype 
qs = da.findGenotypeLabel('WT', 'sqh') ;
% Now we can see what data lies in this slice via
qs.meta
% and we can extract that data via 
qs.getData()
% The data is stored in a property
qsdata = qs.data ;

%% We can grab all the data associated with a given embryo
embryoID = '201905091604' ;
qs = da.findEmbryo(embryoID) ;
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
% To find embryos with a t=10 +/- 2 min
snaps = mapWT.findTime(30, 5) 
% To find Runt stains with a t=10 +/- 2 min
runtsnaps = mapWT.findLabelTime('Runt', 20, 2)

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

%% Create the gradients of pair rule genes for tensor analysis
da.makeGradientImages()
da.makeGradientImages({'Runt'})
% Note: Optionally, declare what scales to examine (in pixels)
% sigmas = [5, 10, 20, 30, 50] ;
% steps = [1] ;
% da.makeGradientImages(sigmas, steps)

%% Use dynamic datasets within the lookupMap to build master timeline
% To control how this is performed, toggle da.timeLineMethod
% Align dynamic runt nanobody data against each other 
Options.apCijFrac = 0.0 ;
Options = struct();
da_1.makeMasterTimeline('WT', 'Runt', Options)
%da.makeMasterTimeline('Tlrm9','Eve',Options)

%% Timestamp other data against the master timeline
% To control how this is performed, toggle da.timeStampMethod
% Here timestamp Runt fixed samples against runt timeline
% You can pass options through a struct if you wish.
Options = struct() ;
da.timeStamp('WT', 'Runt', Options)

% Compare leading edge of Runt to TRAILING edge of Eve 
Options.timelineLabel = 'Eve' ;
Options.timelineLeadingTrailing = 'trailing' ;
da.timeStamp('WT', 'Runt', Options)


