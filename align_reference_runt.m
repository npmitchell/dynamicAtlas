% align_reference_runt.m 
% Script for aligning dynamic runt nanobody data against each other
% 
% NPMitchell 2020

%% ADD PATHS
% Where is the dynamic data?
atlasDir = '/Users/npmitchell/Desktop/tmp/' ;
if ~exist(atlasDir)
    atlasDir = '/mnt/crunch/Atlas_Data/Atlas_Data/' ;
end
% Where are git repos?
gitDir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/' ;
if ~exist(gitDir, 'dir')
    gitDir = '/data/code/' ;
end
gutDir = fullfile(gitDir, 'gut_matlab') ;
basicsDir = fullfile(gutDir, 'basics') ;
tiffDir = fullfile(gutDir, 'tiff_handling') ;
plottingDir = fullfile(gutDir, 'plotting') ;

% Where is the time alignment repo?
tlaDir = '/Users/npmitchell/Box/Flies/code/time_alignment_2020/time_align_embryos';
if ~exist(tlaDir, 'dir')
    tlaDir = '/data/code/time_align_embryos/' ;
end
addpath(fullfile(tlaDir, 'nanconv')) ;
fmDir = fullfile(gutDir, 'toolbox_fast_marching/toolbox_fast_marching/') ;
fmDir2 = fullfile(fmDir, 'mex') ;
fmDir3 = fullfile(fmDir, 'toolbox') ;
addpath(basicsDir) ;
addpath(tiffDir) ;
addpath(plottingDir) ;
addpath(tlaDir) ;
addpath(fmDir) ;
addpath(fmDir2) ;
addpath(fmDir3) ;

%% OPTIONS
% Save each runt nanobody (curated) MIP as ./date/cylinder1_max.tif
runtNBodyDir = fullfile(atlasDir, 'Runt-Nanobody') ;
% where should time alignment and calibration dirs be dumped?
alignDir = atlasDir ;
mipfn = 'cylinder1_max.tif' ;
outdir = fullfile(atlasDir, 'Runt-Nanobody_time_alignment') ;
ssfactor = 4 ;              % subsampling factor before computing corr
% Correlation options
corr_method = 'realspace' ; % realspace or phase method for correlation 
                            % If realspace, does not tranlate, phase allows
                            % dx,dy then computes realspace corr on shifted
                            % image.
stripe7corr_method = 'dist' ; 
hard = 4 ;                  % which experiment is the master timeline -- hard is its index
corrOutDir = fullfile(outdir, sprintf([corr_method '_corr_ss%02d'], ssfactor)) ;

dirs2make = {outdir, corrOutDir} ;
for ii = 1:length(dirs2make)
    dir2make = dirs2make{ii} ;
    if ~exist(dir2make, 'dir')
        mkdir(dir2make)
    end
end

% General options
preview = false ;           % display intermediate results
overwrite = false ;         % overwrite previous results
thres = 0.5 ;

% define naming
mipfnBase = mipfn(1:end-4) ;
extn = [sprintf('_ss%02d', ssfactor) '_' corr_method] ;

exptnames = {'202001141730', '202001141943', '202001142033', ...
    '202001150004', '202001210000', '202001210044'} ;
expts = cell(length(exptnames), 1) ;
for ii = 1:length(exptnames)
    expts{ii} = fullfile(runtNBodyDir, exptnames{ii}) ;
end
% expts = subdirs(runtNBodyDir) ;
substr = sprintf('_ss%02d', ssfactor) ;
ssmipfn = [mipfnBase substr '.tif'];

%% Aesthetics
colorset = define_colors() ;
blue = colorset(1, :) ;
red = colorset(2, :) ;
yellow = colorset(3, :) ;
purple = colorset(4, :) ;
green = colorset(5, :) ;
sky = colorset(6, :) ;
colors = [yellow; sky] ;
gray = [0.5, 0.5, 0.5] ;

%%
align_reference_timing
