function [DSCdata, sortidx] = DSCtools_sortByMass(DSCdata, mode)
% DSCdata = DSCtools_sortByMass(DSCdata)
%
% Sorts DSCdata struct array by sample mass (field "mass")
%
% INPUT:   DSCdata --> DSC data structure as returned by DSCtools_readFile
%             mode --> 'ascend' or 'descend', passed to sort()   [default: 'ascend']
%
% OUTPUT:  DSCdata --> sorted DSC data structure
%          sortidx --> indices as returned by sort();
%
% Author: Andreas Sommer, Jun2026
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu

% mode specified?
if (nargin < 2)
   mode = 'ascend';
end

% extract the rates and sort
masses = [DSCdata.mass];
[~, sortidx] = sort(masses, mode);
DSCdata = DSCdata(sortidx);

end

  