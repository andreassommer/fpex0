function string = DSC204_replaceDecimalDelimiter(string, enc)
   % string = DSC204_replaceDecimalDelimiter(string, enc)
   %
   % Replaces the decimal delimiter as specified in enc.
   % If no encoding is given, the default is used.
   %
   % INPUT:   string --> string or cell array of string to process
   %             enc --> encoding structure
   % 
   % OUTPUT:  string --> processed string(s)
   %
   % Andreas Sommer, Mar2017
   % andreas.sommer@iwr.uni-heidelberg.de
   % email@andreas-sommer.eu
   %
   
   % check if enc is given
   if (nargin<2)
      enc = DSC204_getEncodingDefaults();
   end
   
   % replace the decimal delimiter
   string  = strrep(string, enc.decimalDelimiter, '.');
   
end