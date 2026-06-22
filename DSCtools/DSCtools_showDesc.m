function DSCtools_showDesc(DSCdata, field)
% DSCtools_showDesc(DSCdata, field)
%
% Displays the description of the DSC data structures stored in DSCarray.
% 
% INPUT:  DSCdata --> DSC data structure as returned by DSCtools_readFile()
%           field --> field name to be shown
%
% OUTPUT: on display
%
% Author: Andreas Sommer, Jun2026
% email@andreas-sommer.eu

% get list of all filenames
filenames   = arrayfun(@(x) x.rawData.fileSpec, DSCdata, 'UniformOutput', false);
filelengths = cellfun(@length, filenames);
maxfilelen  = max(filelengths);

% get list of all IDs
maxidlen = max( cellfun(@length, {DSCdata.ID}) );

% walk through DSC data and display
for i = 1:length(DSCdata)
   di = DSCdata(i);
   filespec = di.rawData.fileSpec;
   id       = di.ID;
   desc     = di.desc;
   fprintf('%s\n', getEntry(i, id, maxidlen, desc, field, filespec, maxfilelen));
end

% finito
end


%% HELPERS
function msg = getEntry(i, id, maxidlen, desc, field, file, maxfilelen)
   if isfield(desc, field)
      fieldentry = string(desc.(field));
   else
      fieldentry = '[N/A]';
   end
   msg = sprintf('#%-4d | ID: %*s  |  FILE: %*s  |  field %s: %s', i, maxidlen, id, maxfilelen, file, field, fieldentry);
end
