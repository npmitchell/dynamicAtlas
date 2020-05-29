function [fullPaths, subNames] = subdirs(directory, subdirs2ignore)
% SUBDIRS(directory, subdirs2ignore)
%
%
% NPMitchell 2020

% Unpack options
% choice: always ignore default hidden directories
if nargin > 1
    subdirs2ig = {'.', '..'} ;
    for qq = 1:length(subdirs2ignore)
        subdirs2ig{length(subdirs2ig)+1} = subdirs2ignore{qq} ;
    end
    subdirs2ignore = subdirs2ig ;
else
    subdirs2ignore = {'.', '..'} ;
end

items = dir(directory) ;
dirFlags = [items.isdir] ;
% Extract only those that are directories.
out = items(dirFlags) ;

% Prepare for output as cell array(s)
dmyk = 1 ;
fullPaths = {} ;
if nargout > 1
    subNames = {} ;
end
% Populate cell arrays
for ii = 1:length(out)
    % Ingore hidden paths and other optional ignoreable subdirs
    % To do so, ask if this out(ii).name is in subdirs2ignore
    is_ok = true ;
    for qq=1:length(subdirs2ignore)
        if is_ok && strcmp(out(ii).name, subdirs2ignore{qq}) 
            is_ok = false ;
        end
    end
    
    % add to the output cell if this name is not in the ingored cell array
    if is_ok
        if nargout > 1
            subNames{dmyk} = out(ii).name ;
        end    
        fullPaths{dmyk} = fullfile(out(ii).folder, out(ii).name) ;
        dmyk = dmyk + 1 ;
    end
end