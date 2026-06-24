function DSCtools_showDesc(DSCdata, field)
% DSCtools_showDesc(DSCdata, field)
%
% Displays the description of the DSC data structures stored in DSCarray.
% 
% INPUT:  DSCdata --> DSC data structure as returned by DSCtools_readFile()
%           field --> field name to be shown
%                     if field is empty, a list of all fields is shown.
%
% OUTPUT: on display
%
% Author: Andreas Sommer, Jun2026
% email@andreas-sommer.eu

% special case: field is empty
if isempty(field)
   fieldlist = arrayfun(@(x) fieldnames(x.desc), DSCdata, 'UniformOutput', false);
   fieldlist = vertcat(fieldlist{:});  % collect all
   mfieldlen = length(fieldlist) / length(DSCdata);
   fields    = unique(fieldlist);      % get unique list
   makeMessage('Found: %d unique fields in DSCdata.desc, mean field count is %g.\n', length(fields), mfieldlen);
   for i = 1:length(fields)
      makeMessage('#%4d  %s\n', i, fields{i});
   end
   return
end


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
   makeMessage('%s\n', createOutput(i, id, maxidlen, desc, field, filespec, maxfilelen));
end

% finito
end


%% HELPERS
function msg = createOutput(i, id, maxidlen, desc, field, file, maxfilelen)
   if isfield(desc, field)
      fieldentry = string(desc.(field));
   else
      fieldentry = '[N/A]';
   end
   msg = sprintf('#%-4d | ID: %*s  |  FILE: %*s  |  field %s: %s', i, maxidlen, id, maxfilelen, file, field, fieldentry);
end
