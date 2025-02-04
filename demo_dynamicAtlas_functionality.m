%% Demo DynamicAtlas Functionality
%
% Demo script for the core functions of the Matlab-based atlas code:
% creating an altas structure, creating a master timeline of live data,
% and timestamping fixed data onto the master timeline.
%
% Walkthrough is intended for use on the demo dataset with file name 
% 'DEMO_DATASET.tar.lz4', on the Zenodo repository at the following URL: 
% 
% https://doi.org/10.5281/zenodo.14792464
%
% See corresponding tutorial document included with the publication.
%
% NPMitchell 2020,
% Edited by Vishank Jain-Sharma 2025

%% I. Clear our environment

%clear workspace
clear      

%clear command window
clc      

%close open figures
close all    


%% II. Add paths to atlas data and code to MATLAB search path

%Path to the folder containing the data

%atlasPath = '/PATH/TO/ATLAS/DATA/FOLDER';
%e.g.
atlasPath = '/Users/Vishank/Documents/DynamicAtlas_Demo';

%Path to the folder containing the code

%codePath = '/PATH/TO/ATLAS/CODE/FOLDER';
%e.g.
codePath = '/Users/Vishank/Documents/Dynamic_Atlas_Testing/dynamicAtlas';

%adds atlas path and subfolders to MATLAB search path
addpath(genpath(codePath));

%adds code path and subfolders to MATLAB search path
addpath(genpath(codePath));
pkgDir = fullfile(codePath, '+dynamicAtlas');
cd(pkgDir);
addpath(genpath('./'));
cd(atlasPath)

%% III. Build the dynamicAtlas object
%
% Build dynamic atlas with all genotypes in the atlasPath using:
% da = dynamicAtlas.dynamicAtlas(atlasPath);
%
% Or, choose which genotypes to include in atlas as below.
% By default, all are included

%options specifying how to construct the atlas
Options = struct();       

%method of timeline construction
Options.timeLineMethod = 'realspace';  

%constructing the atlas with wildtype genotypes (WT)
da = dynamicAtlas.dynamicAtlas(atlasPath, {'WT'}, Options); 

%% IV. List the properties of the dynamic atlas

properties(da)

%% V. List the methods of the dynamic atlas

methods(da)

%% VI. Grab all the metadata from a label, and store in a 'queriedSample' 

%data from the wildtype (WT) genotype containing the Runt label
qs = da.findGenotypeLabel('WT','Runt'); 

%% VII. Load in the genotype's data into the queriedSample object

%displays metadata of the queriedSample
qs.meta

%loads in the corresponding data
qs.getData()

%displays the data that was loaded
qs.data

%% VIII. Compute particle-image velocimetry (PIV) on the above data

%options specifying how to compute the PIV
Options = struct();

%overwrite existing PIV computations
Options.overwrite = false;

%computes PIV on the queriedSample data
qs.ensurePIV(Options);

%% IX. Grab all the data associated with a single embryo

%ID of the specific embryo dataset we are querying
embryoID = '202001141730';

%queriedSample with just this embryo's metadata
qs2 = da.findEmbryo(embryoID) ;

%gets the data of this embryo
qs2.getData();

%% X. Use dynamic datasets to build a master timeline
%
% Aligning dynamic run nanobody data against each other

%options specifying how to compute the timeline
Options = struct();

%saving the plots generated whle computing the timeline
Options.save_images = 1;

%makes the master timeline from the live WT datasets with the Runt label
da.makeMasterTimeline('WT','Runt', Options)


%% XI. Timestamp fixed data against the master timeline
%
% Timestamping Runt fixed samples against the master timeline
% Options can be passed through a struct if desired.

%path to the folder of the embryo chosen as the master timeline designee
%(done earlier in the timeline creation block X, the folder will have 
% 'master_timeline_designee.txt' within it)
masterDesigneeDir = '/Users/Vishank/Documents/DynamicAtlas_Demo/WT/Runt/202001141730';

%width(s) of Gaussian(s) to use while smoothing the images
sigmas = [20];
%steps used by the gradient during the computation
steps = [1];
%specifying to compute gradients only on the fixed datasets
%(only fixed are timestamped, so only necessary to compute these here)
fixedOnly = 1;

%makes gradient images of the fixed data to use in the alignment
makeGradientImages(da,'Runt',sigmas,steps,fixedOnly);

%options specifying how to timestamp the fixed samples
Options = struct();
%passes in the directory of the master timeline designee
Options.masterDesigneeDir = masterDesigneeDir;
%indicating that stripe information should be loaded in from the .mat
%variable stored in the folder
Options.loadStripeMat = 1;

%timestamps fixed samples with the Runt label to the master timeline
da.timeStamp('WT', 'Runt', Options) 

disp('Demo done.')