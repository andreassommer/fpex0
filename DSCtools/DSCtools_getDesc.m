function content = DSCtools_getDesc(DSCdata, field)
% content = DSCtools_getDesc(DSCdata, field)
%
% Retrieves fields from the description of the DSC data structures stored in DSCarray.
% 
% INPUT:  DSCdata --> DSC data structure as returned by DSCtools_readFile()
%           field --> field name to be shown
%                     if field is empty, a list of all fields is returned
%
% OUTPUT:  content of field as cellstr (cell array of char arrays).
%
% Author: Andreas Sommer, Jun2026
% email@andreas-sommer.eu

% special case: field is empty
if isempty(field)
   fieldlist = arrayfun(@(x) fieldnames(x.desc), DSCdata, 'UniformOutput', false);
   fieldlist = vertcat(fieldlist{:});  % collect all
   fields    = unique(fieldlist);      % get unique list
   %makeMessage('Found: %d unique fields in DSCdata.desc, mean field count is %g.\n', length(fields), length(fieldlist) / length(DSCdata));
   content = fields;
   return
end

% walk through DSC data and extract -- we don't use arrayfun as a field might be missing in some entries
content = cell(size(DSCdata));
for i = 1:numel(DSCdata)
   desc = DSCdata(i).desc;
   if isfield(desc, field)
      content{i} = desc.(field);
   end
end

% finito
end

