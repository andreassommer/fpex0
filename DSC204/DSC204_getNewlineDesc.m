function [newlinestring, len] = DSC204_getNewlineDesc(fid)
   % [newlinestring, len] = DSC204_getNewlineDesc(fid)
   %
   % Extracts the newline string and its length from open file.
   %
   % Author: Andreas Sommer, Mar2017
   % andreas.sommer@iwr.uni-heidelberg.de
   % email@andreas-sommer.eu
   
   % detect line ending markers; retrieve their length
   lineW = fgets(fid);            % read line with newline chars
   fseek(fid, -length(lineW), 0); % rewind
   line  = fgetl(fid);            % read line without newline chars
   fseek(fid, -length(lineW), 0); % rewind (with newline!)
   
   % set newline string and its length
   len = length(lineW) - length(line);
   newlinestring = lineW(end-len+1:end);

end