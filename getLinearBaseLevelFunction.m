function fun = getLinearBaseLevelFunction(T,M,varargin)
% fun = getLinearBaseLevelFunction(T,M,varargin)
%
% Returns a function delivering the linearly interpolated base-level at a
% certain temperature/time
%
% INPUT:   T:  x-coordinates (temperatures, times)
%          M:  y-coordinates (1D profile of Cp)
%
% Additional configuration can be done via variable input arguments.
%
% OUTPUT:  f(t): delivering linearly interpolated base-level at t
%
% Assumptions:
%   * base level is linear in T
%
% Algorithm: 
%   * leftval = calculate mean value of first-n values starting from value closest to first-T
%   * rightval = calculate mean value of last-n values up to value closest to last-T
%   * make a function that interpolates linearly between leftval and rightval
%
%
% Copyright 2016-2022, Andreas Sommer  code@andreas-sommer.eu
%
% Copying and distribution of this file, with or without modification, are permitted in any medium without royalty, 
% provided the copyright notice and this notice are preserved. This file is offered as-is, without any warranty.



% defaults
firstN = 10;    firstT = T(1);      firstTidx = [];
lastN  = 10;    lastT  = T(end);    lastTidx  = [];

% config 
if olHasOption(varargin,'firstN'),       firstN = olGetOption(varargin,'firstN');    end
if olHasOption(varargin,'firstT'),       firstT = olGetOption(varargin,'firstT');    end
if olHasOption(varargin,'firstTidx'), firstTidx = olGetOption(varargin,'firstTidx'); end
if olHasOption(varargin,'lastN'),         lastN = olGetOption(varargin,'lastN');     end
if olHasOption(varargin,'lastT'),         lastT = olGetOption(varargin,'lastT');     end
if olHasOption(varargin,'lastTidx'),   lastTidx = olGetOption(varargin,'lastTidx');  end

% find indices of first-T and last-T
if isempty(firstTidx), [~,firstTidx] = min(abs(T-firstT)); end
if isempty( lastTidx), [~, lastTidx] = min(abs(T- lastT)); end
   
% calculate mean values
try
   T0base = mean(M(firstTidx+firstN));
   TFbase = mean(M(lastTidx-lastN));
catch err  % probably Index-Out-Of-Bounds
   rethrow(err)
end

% interpolation function
function f = interpfun(t)
   f = interp1([firstT lastT],[T0base TFbase],t,'linear');
end

% return handle to the interpolation function
fun = @interpfun;

% finito
return;

end



   





   
