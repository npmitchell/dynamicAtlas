%% MASTER TIMING PIPELINE
% Template pipeline for loading all pair rule gene patterns
%
% Pair rule genes:
% even-skipped (eve)
% hairy
% odd-skipped (odd)
% paired (prd)
% Runt (runt)
% fushi-Tarazu (ftz)
% odd-paired ()
% sloppy paired (slp) --> also gap gene
% Tenascin major
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
tlaDir = '/Volumes/minimalData/code/';
% tlaDir = '/Users/mattlefebvre/Desktop/Code/code';
% tlaDir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/dynamicAtlas/code';
cd(fullfile(tlaDir)) ;
addpath(genpath('dynamicAtlasCode')) ;
cd('dynamicAtlasCode')
addpath(genpath('+dynamicAtlas'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define atlasPath to be where the dynamicAtlas resides (the parent
% dynamicAtlas directory, not the project directory '+dynamicAtlas')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
atlasPath = '/Volumes/Elements/Atlas_Data' ;

%% Build the dynamicAtlas

% Or choose which labels to include in atlas. Use only WT genotypes
% Note that Even_Skipped and Even_Skipped-YFP are both labels of Eve, but one has a live reporter
pairRuleGenes = {'Runt', 'Even_Skipped', 'Even_Skipped-YFP', 'Hairy', ...
    'Paired', 'Fushi_Tarazu', 'Sloppy_Paired'} ;
options = struct() ;
options.labels = pairRuleGenes ;
da = dynamicAtlas.dynamicAtlas(atlasPath, {'WT'}, options) ;

% Build videos for each PR gene
timestamps = 10:40 ;    % timestamps to seek in each gene
time_unc = 4 ;          % when seeking a timestamp t, seek any in (t-time_unc, t+time_unc).
normalizeEach = true ;  % let each image contribute equally despite brightness variations across samples
preview = false ;       % view intermediate results 

% preallocate cell array for the movies
prMovies = cell(length(pairRuleGenes), 1) ;
for prID = 1:length(pairRuleGenes)
    for tID = 1:length(timestamps)
        qs = findGenotypeLabelTime(da, 'WT', pairRuleGenes{prID}, time, time_unc) ;
        meanIm = qs.getMeanData(normalizeEach, preview) ;
        
        % if this is the first one, instantiate/preallocate the movie array
        if tID == 1
            prMovies{prID} = nan(length(timestamps), ...
                size(meanIm, 1), size(meanIm, 2)) ;
        end
        
        prMovies{prID}(tID, :, :) = meanIm ;
    end
end

% View results
for prID = 1:length(pairRuleGenes)
end

