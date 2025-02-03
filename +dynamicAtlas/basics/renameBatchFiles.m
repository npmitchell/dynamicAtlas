function renameBatchFiles(directory, fnSearchStr, oldchar, newchar, dryrun)
% change (batch) file names
% 
% Parameters
% ----------
% directory : str 
%   path to file to rename
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
% NPMitchell 2020

if nargin < 1
    directory= './';
end
if nargin < 2
    fnSearchStr = '*stripe7*.png' ;
end
if nargin < 3
    oldchar = '.png' ;
end
if nargin < 4
    newchar = '_new.png' ;
end
if nargin < 5
    dryrun = false ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
existingStrReplace= true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if existingStrReplace
    %%%%%%%%%%%%%%%%%%%
    %Replace characters at position of preexisting string in name
    toRename = dir(fullfile(directory, fnSearchStr)) ;
    if ~isempty(toRename)
        nchar = strfind(toRename(1).name, oldchar)  ;
        if length(nchar) > 1
            nchar = nchar(1) ;
        end

        lenstr= length(toRename(1).name);
        fprintf('First new name is:\n  %s\n',[toRename(1).name(1:nchar-1), newchar, toRename(1).name(nchar+length(oldchar):lenstr)])
            
        for ii=1:length(toRename)
            % update the character index
            nchar = strfind(toRename(ii).name, oldchar) ;
            if size(nchar) > 1
                nchar = nchar(1) ;
                replace_file = true ;
            elseif size(nchar) == 0
                replace_file = false ;
            else
                replace_file = true ;
            end

            if replace_file
                lenstr= length(toRename(ii).name);
                newname = [toRename(ii).name(1:nchar-1), newchar,...
                    toRename(ii).name(nchar+length(oldchar):lenstr)] ;
                if dryrun
                    disp(fullfile(directory, toRename(ii).name)) 
                    disp( [ ' > ' fullfile(directory, newname)])
                else
                    disp( [ ' source > ' fullfile(directory, toRename(ii).name)]) 
                    disp( [ '   dest > ' fullfile(directory, newname)])
                    movefile(fullfile(directory, toRename(ii).name), ...
                        fullfile(directory, newname))
                end
            end
        end
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Replace characters at known position in string
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% knownLocReplace = false ;
% pos = 10 ;  % position of characters after which to replace
% if knownLocReplace
%     %Replace characters at known position in string
%     toRename = dir(fnSearchStr);
%     nchar = length(toRename(1).name) ;
%     fprintf('First new name is:\n  %s\n',[toRename(ii).name(1:pos), newchar])
% 
%     for ii=1:length(toRename)
%         lenstr= length(toRename(ii).name);
%         newname = [toRename(ii).name(1:pos), newchar] ;
%         movefile([directory, toRename(ii).name], [directory, newname])
%     end
% end

