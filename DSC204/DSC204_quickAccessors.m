function [t, T, uV, sf, mW] = DSC204_quickAccessors(dataMat, Tmin, Tmax)
% [t, T, uV, sf, mW] = DSC204_quickAccessors(dataMat, Tmin, Tmax);
%
% Separates the dataMat array of an DSC204data structure as returned by DSC204_readFile.
%
% INPUT:  DSC204data --> dataMat array
%               Tmin --> minimum temperature to restrict data to   (default: -inf)
%               Tmax --> maximum temperature to restrict data to   (default: +inf)
%
% OUTPUT:  t --> times
%          T --> temperatures (of reference)
%         uV --> DSC microvolt signal
%         sf --> sensitivity factors
%         mW --> heat flux in milli-Watts
%
% Author:  Andreas Sommer, Mar2017, Aug2022
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu

% default indices
idx_t  = 2; 
idx_T  = 1;
idx_uV = 3;
idx_sf = 4;

% import settings
global DSC204settings
idx_t  = getSetting(DSC204settings, 'dataMAT_column_t' , idx_t);
idx_T  = getSetting(DSC204settings, 'dataMAT_column_T' , idx_T);
idx_uV = getSetting(DSC204settings, 'dataMAT_column_uV', idx_uV);
idx_sf = getSetting(DSC204settings, 'dataMAT_column_sf', idx_sf);

% check input arguments
if (nargin < 3); Tmax = +inf; end
if (nargin < 2); Tmin = -inf; end

% quick accessors
t  = dataMat(:,idx_t);
T  = dataMat(:,idx_T);
uV = dataMat(:,idx_uV);
sf = dataMat(:,idx_sf);

% restrict to specified range
if ~all(isinf([Tmin Tmax]))
   idx = ((T >= Tmin) & (T <= Tmax));
   t  = t(idx);
   T  = T(idx);
   uV = uV(idx);
   sf = sf(idx);
end

% build heat flux if requested
if (nargout >= 5)
   mW = uV ./ sf;
end

   
end