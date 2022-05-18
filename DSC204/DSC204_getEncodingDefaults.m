function enc = DSC204_getEncodingDefaults()
   % defaults = DSC204_getEncodingDefaults()
   %
   % Defaults like encodings, delimiters, etc. are stored here.
   %
   % Author: Andreas Sommer, Mar2017
   % andreas.sommer@iwr.uni-heidelberg.de
   % email@andreas-sommer.eu
   %
   
   % defaults
   enc.machineformat    = 'native';
   enc.encoding         = 'ISO-8859-1';   % this one fits best, but still does not get "degree" and "micro"
   enc.fieldDelimiter   = ';';
   enc.decimalDelimiter = ',';

   
   % depending on the current character set, use the following substitutions
   try
      characterSet = feature('DefaultCharacterSet');
   catch
      characterSet = 'unknown';
   end
   unicode_deg = native2unicode([176],'ISO-8859-1');       %#ok<NBRAK>
   unicode_mu  = native2unicode([181],'ISO-8859-1');       %#ok<NBRAK>
   switch upper(characterSet)
      case 'UTF-8'
         char_deg = native2unicode(unicode2native(unicode_deg), 'UTF-8');
         char_mu  = native2unicode(unicode2native(unicode_mu) , 'UTF-8');
      case 'US-ASCII'
         char_deg = 'deg';    % does not exist in 7-bit US-ASCII
         char_mu  = native2unicode(unicode2native(unicode_mu) , 'US-ASCII');
      case 'WINDOWS-1252'
         char_deg = native2unicode(unicode2native(unicode_deg), 'WINDOWS-1252');
         char_mu  = native2unicode(unicode2native(unicode_mu) , 'WINDOWS-1252');
      case 'ISO-8859-1'
         char_deg = native2unicode(unicode2native(unicode_deg), 'ISO-8859-1');
         char_mu  = native2unicode(unicode2native(unicode_mu) , 'ISO-8859-1');
      otherwise
         char_deg = 'deg';
         char_mu  = 'micro';
   end
   % substitutions    = {'deg'    , 'micro'  };
   enc.subst_codes    = { 155     ,   145    };  % these are the codes found in the file
   enc.subst_strings  = { char_deg,  char_mu };  % must be adjusted to encoding

   % store these characters
   enc.degC = char_deg;
   enc.mu   = char_mu;
  
end