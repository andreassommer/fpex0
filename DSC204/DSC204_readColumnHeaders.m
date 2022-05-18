function columnHeaders = DSC204_readColumnHeaders(fid, enc)
   % columnHeaders = DSC204_readColumnHeaders(fid, enc)
   %
   % Searches and returns the column headers in an DSC204 export file.
   %
   % INPUT:   fid --> file id to DSC204 file
   %          enc --> encoding structure (defaults)
   %
   % OUTPUT:  columnHeaders --> cell array of strings
   %
   % Author: Andreas Sommer, Mar2017
   % andreas.sommer@iwr.uni-heidelberg.de
   % email@andreas-sommer.eu
   %
   
   
   % save current file position, rewind file
   filePosition = ftell(fid);
   frewind(fid);
   
   % prepare empty headers
   columnHeaders = {};
   
   % walk through file and search column header tag ##
   while ~feof(fid)
      % read line (with newline chars)
      line = fgets(fid);
      linelen = length(line);
      % check if column header tag ## is found
      if (linelen>=2) && strcmp(line(1:2),'##')
         line = DSC204_replaceChars(line, enc);  % replace strange characters in file header
         columnHeaders = textscan(line(3:end),'%s','Delimiter',enc.fieldDelimiter);
         columnHeaders = columnHeaders{1}; % unwrap first cell layer
         break
      end
   end
   
   % check if eof
   if feof(fid), warning('Unexpected end of file!'); end
   
   % ensure column headers form a single row
   columnHeaders = reshape(columnHeaders,1,[]);
   
   % finito: restore file position
   fseek(fid, filePosition, 'bof');
   return

   
   end