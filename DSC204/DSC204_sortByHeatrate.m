function [DSCdata, sortidx] = DSC204_sortByHeatrate(DSCdata, mode)
% DSCdata = sortByRate(DSCdata)
%
% Sorts DSCdata struct array by heatrate (field "rate")
%
% INPUT:   DSCdata --> DSC data structure as returned by DSC204_readFile
%             mode --> 'ascend' or 'descend', passed to sort()
%
% OUTPUT:  DSCdata --> sorted DSC data structure
%          sortidx --> indices as returned by sort();
%
% Author: Andreas Sommer, Aug2022
% andreas.sommer@iwr.uni-heidelberg.de
% email@andreas-sommer.eu

% mode specified?
if (nargin < 2)
   mode = 'ascend';
end

% extract the rates and sort
rates = [DSCdata.rate];
[~, sortidx] = sort(rates, mode);
DSCdata = DSCdata(sortidx);

end

  