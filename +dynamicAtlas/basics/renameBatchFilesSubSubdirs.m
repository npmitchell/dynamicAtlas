function renameBatchFilesSubSubdirs(directory, fnSearchStr, oldchar, newchar, dryrun)
% change (batch) file names in every subdir of a directory
% 
% Parameters
% ----------
% directory : str 
%   path of grandparent directory of file to rename
% fnSearchStr : str
%   string specifier to find files to rename
% oldchar : str
%   string to replace in old filenames
% newchar : str 
%   replacement string in new filenames
% dryrun : optional bool, default=false
%   do not actually replace files, just test I/O
%
% Returns 
% -------
% <none>
%
% Example Usage
% -------------
% direc = '/Volumes/minimalData/Atlas_Data/WT/' ;
% fnSearchStr = 'timematch_curve7_*.*' ;
% oldchar = 'timematch_curve7_' ;
% newchar = 'timematch_Runt_Runtstripe7_' ;
% dryrun = true ;
% renameBatchFilesSubSubdirs(direc, fnSearchStr, oldchar, newchar, dryrun)
%
% NPMitchell 2020

subdir2do = subdirs(directory) ;
for qq = 1:length(subdir2do)
    subdir = subdir2do{qq} ;
    renameBatchFilesSubdirs(subdir, fnSearchStr, oldchar, newchar, dryrun)
end