function string = DSC204_replaceChars(string, enc)
   % string = DSC204_replaceChars(string, enc)
   %
   % Replaces characters in string as specified in encoding structure enc.
   %
   % INPUT:   string --> string to process
   %             enc --> encoding structure
   %
   % OUTPUT:  string --> processed string
   %
   % Andreas Sommer, Mar2017
   % andreas.sommer@iwr.uni-heidelberg.de
   % email@andreas-sommer.eu
   %
   
   for k = 1:length(enc.subst_codes)
      string = strrep(string, char(enc.subst_codes{k}), enc.subst_strings{k});
   end
   
end
