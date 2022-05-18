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
if hasOption(varargin,'firstN'),       firstN = getOption(varargin,'firstN');    end
if hasOption(varargin,'firstT'),       firstT = getOption(varargin,'firstT');    end
if hasOption(varargin,'firstTidx'), firstTidx = getOption(varargin,'firstTidx'); end
if hasOption(varargin,'lastN'),         lastN = getOption(varargin,'lastN');     end
if hasOption(varargin,'lastT'),         lastT = getOption(varargin,'lastT');     end
if hasOption(varargin,'lastTidx'),   lastTidx = getOption(varargin,'lastTidx');  end

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



   





   
