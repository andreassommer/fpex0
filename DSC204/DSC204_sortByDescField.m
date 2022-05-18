function [DSCdata, sortidx] = DSC204_sortByDescField(DSCdata, fieldname, transformation)
% DSCdata = sortByTstep(DSCdata, fieldname)
%
% Sorts DSCdata struct array by .desc.(fieldname) field.
%
% INPUT:   DSCdata --> DSC data structure as returned by DSC204_readFile
%            field --> field name to sort after
%   transformation --> function handle that will be applied upon each data
%                      string in specified field before the sorting (e.g. @str2num)
%                      leave empty for none
%
% OUTPUT:  DSCdata --> sorted DSC data structure
%          sortidx --> indices as returned by sort();
%
%
% Examples:  dsc2 = DSC204_sortByDescField(dsc,'SAMPLEMASSmg',@(x) str2num(DSC204_replaceDecimalDelimiter(x,dsc2(1).fileEnc)))
%             
%
% Author: Andreas Sommer, Mar2017
% andreas.sommer@iwr.uni-heidelberg.de
% email@andreas-sommer.eu

% transformation specified?
if (nargin < 3); transformation = []; end

% extract the Tsteps
desc   = {DSCdata.desc};
fields = cellfun(@(x) x.(fieldname), desc, 'UniformOutput', false);

% apply transformation
if ~isempty(transformation)
   fields = cellfun(@(x) transformation(x), fields, 'UniformOutput', false);
   % check if values are numeric, then transform the cell array into a vector
   if all(cellfun(@isnumeric,fields)) && all(cellfun(@isscalar,fields))
      fields = cell2mat(fields);
   end
end

% sort
[~, sortidx] = sort(fields);
DSCdata = DSCdata(sortidx);

end

  