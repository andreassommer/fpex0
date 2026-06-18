function [DSCdata, sortidx] = DSCtools_sortByFileSpec(DSCdata)
% DSCdata = DSCtools_sortByFileSpec(DSCdata)
%
% Sorts DSCdata struct array by .fileSpec field, in lexicographic order.
%
% INPUT:   DSCdata --> DSC data structure as returned by DSCtools_readFile(s)
%
% OUTPUT:  DSCdata --> sorted DSC data structure
%          sortidx --> indices as returned by sort();
%
% Author: Andreas Sommer, Mar2017
% andreas.sommer@iwr.uni-heidelberg.de
% email@andreas-sommer.eu

% extract the Tsteps
filespecs = {DSCdata.fileSpec};

% sort
[~, sortidx] = sort(filespecs);
DSCdata = DSCdata(sortidx);

end

  