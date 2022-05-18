function [DSCdata, sortidx] = DSC204_sortByTstep(DSCdata, mode)
% DSCdata = sortByTstep(DSCdata)
%
% Sorts DSCdata struct array by .Tinfo.Tstep field
%
% INPUT:   DSCdata --> DSC data structure as returned by DSC204_readFile
%             mode --> 'ascend' or 'descend', passed to sort()
%
% OUTPUT:  DSCdata --> sorted DSC data structure
%          sortidx --> indices as returned by sort();
%
% Author: Andreas Sommer, Mar2017
% andreas.sommer@iwr.uni-heidelberg.de
% email@andreas-sommer.eu

% mode specified?
if (nargin < 2)
   mode = 'ascend';
end

% extract the Tsteps
Tinfos = {DSCdata.Tinfo};
Tsteps = cellfun(@(x) x.Tstep, Tinfos);

% sort
[~, sortidx] = sort(Tsteps, mode);
DSCdata = DSCdata(sortidx);

end

  