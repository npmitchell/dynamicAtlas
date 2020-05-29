%Noah Mitchell
%change (batch) file names

%make into a function...

directory= './';
fnSearchStr = '*stripe7*.png' ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
existingStrReplace= true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if existingStrReplace
    %%%%%%%%%%%%%%%%%%%
    %Replace characters at position of preexisting string in name
    oldchar = '.png' ;
    newchar = '_dist.png' ;
    toRename = dir(fullfile(directory, fnSearchStr));
    nchar = strfind(toRename(1).name, oldchar) ;

    lenstr= length(toRename(1).name);
    fprintf('First new name is:\n  %s\n',[toRename(1).name(1:nchar-1), newchar, toRename(1).name(nchar+length(oldchar):lenstr)])

    for ii=1:length(toRename)
        lenstr= length(toRename(ii).name);
        newname = [toRename(ii).name(1:nchar-1), newchar, toRename(ii).name(nchar+length(oldchar):lenstr)];
        movefile([directory, toRename(ii).name], [directory, newname])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Replace characters at known position in string
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
knownLocReplace = false ;

if knownLocReplace
    %Replace characters at known position in string
    newchar = '.tif' ;
    toRename = dir([directory '/20*']);
    nchar = length(toRename(1).name) ;
    fprintf('First new name is:\n  %s\n',[toRename(ii).name(1:nchar-1), newchar])

    for ii=1:length(toRename)
        lenstr= length(toRename(ii).name);
        newname = [toRename(ii).name(1:nchar-1), newchar] ;
        movefile([directory, toRename(ii).name], [directory, newname])
    end
end

