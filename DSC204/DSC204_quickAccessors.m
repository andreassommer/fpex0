function [t, T, uV, sf, mW] = DSC204_quickAccessors(DSC204data, Tmin, Tmax)
% [t, T, uV, sf, mW] = DSC204_quickAccessors(DSC204data, Tmin, Tmax);
%
% Separates the data array of an DSC204data structure as returned by DSC204_readFile.
%
% INPUT:  DSC204data --> structure as returned by DSC204_readFile
%               Tmin --> minimum temperature to restrict data to   (default: -inf)
%               Tmax --> maximum temperature to restrict data to   (default: +inf)
%
% OUTPUT:  t --> times
%          T --> temperatures (of reference)
%         uV --> DSC microvolt signal
%         sf --> sensitivity factors
%         mW --> heat flux in milli-Watts
%
% Author:  Andreas Sommer, Mar2017
% andreas.sommer@iwr.uni-heidelberg.de
% email@andreas-sommer.eu


% check input arguments
if (nargin < 3); Tmax = +inf; end;
if (nargin < 2); Tmin = -inf; end;

% quick accessors
t  = DSC204data.data(:,2);
T  = DSC204data.data(:,1);
uV = DSC204data.data(:,3);
sf = DSC204data.data(:,4);

% restrict to specified range
if ~all(isinf([Tmin Tmax]))
   idx = ((T >= Tmin) & (T <= Tmax));
   t  = t(idx);
   T  = T(idx);
   uV = uV(idx);
   sf = sf(idx);
end

% build heat flux
if (nargout >= 5)
   mW = uV ./ sf;
end

   
end