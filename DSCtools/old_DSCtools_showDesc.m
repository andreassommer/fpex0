function DSCtools_showDesc(DSCarray)
% DSCtools_showDesc(DSCarray)
%
% Displays the description of the DSC-data structures stored in DSCarray.
% 
% INPUT:  DSCarray --> struct array of DSC structures (as returned by DSCtools_readFile or DSCtools_loadall)
%
% OUTPUT: on display
%
% Author: Andreas Sommer, Mar2017
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu



% get list of all description fields 
descfields = arrayfun(@(x) fieldnames(x.desc), DSCarray, 'UniformOutput', false);
descfields = vertcat(descfields{:});       % combine them
descfields = unique(descfields, 'stable'); % keep unique, and no sorting of fields


% now display grouped by field names
for n = 1:length(descfields)
   fieldname = descfields{n};
   makeMessage('\n\nField %s:\n', fieldname);
   % walk through all DSC measurements
   count = 1;
   for k = 1:length(DSCarray)
      try
         fileSpec = DSCarray(k).fileSpec;
         content  = DSCarray(k).desc.(fieldname);
         makeMessage('  %50s: %s\n', fileSpec, content);
         count = count + 1;
      catch err
         makeMessage('Ignoring error @ #%d for file %s: %s\n', n, fileSpec, err.message)
      end
   end
   makeMessage('Count: %d\n', count);
end


% finito
end