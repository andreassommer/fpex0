function dscdata = DSC204_addQuickAccessors(dscdata, Tmin, Tmax)
% dscdata = DSC204_addQuickAccessors(dscdata, Tmin, Tmax)
%
% Separates the raw dataMat array of an DSC204data structure as returned by DSC204_readFile
% into the respective column data.
%
% INPUT:  dscdata --> (vector of) DSC204data structures as returned by DSC204_readFile
%            Tmin --> minimum temperature to restrict data to   (default: -inf)
%            Tmax --> maximum temperature to restrict data to   (default: +inf)
%
% OUTPUT: dscdata --> (vector of) DSC204data with .data field containing:
%            data.t  --> times
%            data.T  --> temperatures (of reference)
%            data.uV --> DSC microvolt signal
%            data.sf --> sensitivity factors
%            data.mW --> heat flux in milliwatts (calculated as uV/sf)
%
% Author:  Andreas Sommer, Mar2017, Aug2022
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu


% check input arguments
if (nargin < 3); Tmax = +inf; end
if (nargin < 2); Tmin = -inf; end

% make function applicable for struct arrays
if length(dscdata) > 1
   dscdata = arrayfun(@(x) DSC204_addQuickAccessors(x, Tmin, Tmax));
   return
end

% extract and store quick accessors
[t, T, uV, sf, mW] = DSC204_quickAccessors(dscdata.rawData.dataMat, Tmin, Tmax);
dscdata.data.t  = t;
dscdata.data.T  = T;
dscdata.data.uV = uV;
dscdata.data.sf = sf;
dscdata.data.mW = mW;

% finito
return
   
end