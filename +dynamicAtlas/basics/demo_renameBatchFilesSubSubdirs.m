%% DEMO FOR renameBatchFilesSubSubdirs

% here, rename all timematch files to include labels
direc = '/Volumes/minimalData/Atlas_Data/WT/' ;
fnSearchStr = 'timematch_curve7_*.*' ;
oldchar = 'timematch_curve7_' ;
newchar = 'timematch_Runt_Runtstripe7_' ;
dryrun = false ; % whether to actually rename files or test I/0
renameBatchFilesSubSubdirs(direc, fnSearchStr, oldchar, newchar, dryrun)

% rename the fit images (pngs) to include the labels too 
dryrun = true ;  % whether to actually rename files or test I/0
fnSearchStr = 'stripe7_chisq_fit.png' ;
oldchar = 'stripe7_' ;
newchar = 'timematch_Runt_Runtstripe7_' ;
renameBatchFilesSubSubdirs(direc, fnSearchStr, oldchar, newchar, dryrun)
