function dscdata = DSCtools_addQuickAccessors(dscdata, Tmin, Tmax)
% dscdata = DSCtools_addQuickAccessors(dscdata, Tmin, Tmax)
%
% Separates the raw dataMat array of an DSCTOOLS_GLOBAL_SETTINGS structure as returned by DSCtools_readFile
% into the respective column data.
%
% INPUT:  dscdata --> (vector of) DSCTOOLS_GLOBAL_SETTINGS structures as returned by DSCtools_readFile
%            Tmin --> minimum temperature to restrict data to   (default: -inf)
%            Tmax --> maximum temperature to restrict data to   (default: +inf)
%
% OUTPUT: dscdata --> (vector of) DSCTOOLS_GLOBAL_SETTINGS with .data field containing:
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
   dscdata = arrayfun(@(x) DSCtools_addQuickAccessors(x, Tmin, Tmax));
   return
end

% extract and store quick accessors -- does not add to memory footprint unless modified!
[t, T, uV, sf, mW] = DSCtools_quickAccessors(dscdata.rawData.dataMat, Tmin, Tmax);
dscdata.data.t  = t;
dscdata.data.T  = T;
dscdata.data.uV = uV;
dscdata.data.sf = sf;
dscdata.data.mW = mW;

% finito
return
   
end