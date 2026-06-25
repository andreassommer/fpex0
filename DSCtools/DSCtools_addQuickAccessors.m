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
% Author:  Andreas Sommer, Mar2017, Aug2022, Jun2026
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

% there should always be t (time) and T (temperature) available
cHeaders = dscdata.rawData.columnHeaders;
[ idx_t  , units_t   ] = DSCtools_identifyColumns(cHeaders, 'Time');  % must be there
[ idx_T  , units_T   ] = DSCtools_identifyColumns(cHeaders, 'Temp');  % must be there
[ idx_DSC, units_DSC ] = DSCtools_identifyColumns(cHeaders, 'DSC');   % must be there
[ idx_sf , units_sf  ] = DSCtools_identifyColumns(cHeaders, 'Sens');  % should be there

% quick access to data
dataMat = dscdata.rawData.dataMat;

% get Temperature - restricted to range
T = dataMat(:, idx_T);
idxTrange = (T >= Tmin) & (T <= Tmax);
T = T(idxTrange);

% get other quantities
t   = dataMat(idxTrange, idx_t);    % time
DSC = dataMat(idxTrange, idx_DSC);  % actual measurement
sf  = dataMat(idxTrange, idx_sf);   % sensitivity [might be empty]

% if sensitivity is empty, set it all to zero
if isempty(sf)
   sf = zeros(size(T));
end

% We now have to interpret the unit of "DSC" -- might be uV, uV/mg, mW, mW/mg
DSC_contains_mW = contains(units_DSC, 'mW');
DSC_contains_uV = contains(units_DSC, 'uV');
if ~xor(DSC_contains_uV, DSC_contains_mW)
   error('Found both units "mW" and "uV" or none of them in DSC. Abort.')
end
if DSC_contains_mW
   mW = DSC;
   uV = mW ./ sf;
end
if DSC_contains_uV
   uV = DSC;
   mW = uV .* sf;
end

% check the other units
DSCtools_checkUnits('Time'       , units_t);
DSCtools_checkUnits('Temperature', units_T);
DSCtools_checkUnits('Scaling'    , units_sf);


% store quick accessors
dscdata.data.t  = t;
dscdata.data.T  = T;
dscdata.data.uV = uV;
dscdata.data.sf = sf;
dscdata.data.mW = mW;

% finito
return
   
end