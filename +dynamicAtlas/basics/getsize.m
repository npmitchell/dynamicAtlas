function getsize(this) 
   %GETSIZE Show the size on disk of an object.
   %   Only works if the properties GetAccess method is set to public
   %
   % NPM 2019
   
   % first try the simple thing
   whos(this)
   % now try the fancy thing
   props = properties(this); 
   totSize = 0; 
   
   for ii=1:length(props) 
      currentProperty = getfield(this, char(props(ii))); 
      s = whos('currentProperty'); 
      totSize = totSize + s.bytes; 
   end
  
   fprintf(1, '%d bytes\n', totSize); 
end