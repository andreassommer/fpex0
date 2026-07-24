function [DSCdata, sortidx] = DSCtools_sortByFileSpec(DSCdata, mode)
% DSCdata = DSCtools_sortByFileSpec(DSCdata)
%
% Sorts DSCdata struct array by .fileSpec field, in lexicographic order.
%
% INPUT:   DSCdata --> DSC data structure as returned by DSCtools_readFile(s)
%             mode --> 'ascend' or 'descend'      [default: 'ascend']
%
% OUTPUT:  DSCdata --> sorted DSC data structure
%          sortidx --> indices as returned by sort();
%
% Author: Andreas Sommer, Mar2017
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu

% sortmode -- sort can only sort cell array ascending, so we have to cope here
if (nargin < 2)
   sortmode = 'a';  % ascend
else
   mode = convertStringsToChars(mode);
   sortmode = mode(1);
   if ~ismember(sortmode, 'ad')
      error('Unknown sort mode: "%s". Use lower case "ascend" or "descend".', mode)
   end
end

% extract the file specs
filespecs = arrayfun(@(x) x.rawData.fileSpec, DSCdata, 'UniformOutput', false);

% sort
[~, sortidx] = sort(filespecs);  % sort can sort cells only ascending!
if (sortmode == 'd')
   sortidx = flip(sortidx);      % reverse sorting if sortmode is descending
end
DSCdata = DSCdata(sortidx);

end

  