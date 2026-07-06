function [allOK, tflist] = DSCtools_assertTrange(DSCsample, Tmin, Tmax, varargin)
% [allOK, tflist] = DSCtools_assertTrange(DSCsampleT, Tmin, Tmax, key-value-pairs)
%
% Checks if specified temperature range is fully available inside DSCsample
%
% INPUT:  DSCsample --> DSC data structure as returned by DSCtools_readFile
%              Tmin --> minimum temperature to restrict data to   [default: -inf]
%              Tmax --> maximum temperature to restrict data to   [default: +inf]
%
% OUTPUT:     allOK --> true if [Tmin, Tmax] is a subset of all available measurements
%            tflist --> tflist(i) is true if [Tmin, Tmax] is a subset of measurements in DSCsample(i),
%                       otherwise false
%
% Author: Andreas Sommer, Jun2026
% code@andreas-sommer.eu

% default args
if (nargin < 3), error('Not enough arguments'); end

% check optional arguments
args = varargin;
[ warnOnFail, args] = olGetOption(args,  'warnOnFail', true );
[errorOnFail, args] = olGetOption(args, 'errorOnFail', false);
olWarnIfNotEmpty(args);

% check for every element of DSCsample, if [Tmin,Tmax] is a subset
if isinf(Tmin), tflist_Tmin = true(size(DSCsample));
else
   tflist_Tmin = arrayfun( @(x) (Tmin >= x.data.T( 1 )) , DSCsample);
end
if isinf(Tmax)
   tflist_Tmax = true(size(DSCsample));
else
   tflist_Tmax = arrayfun( @(x) (Tmax <= x.data.T(end)) , DSCsample);
end
tflist = tflist_Tmin & tflist_Tmax;
allOK = all(tflist);

% process warning and/or error
if ~allOK
   keyboard
   if  warnOnFail, warning( make_failure_message(tflist, Tmin, Tmax) ); end
   if errorOnFail,   error( make_failure_message(tflist, Tmin, Tmax) ); end
end

% finito
return
   
end % of function


%% HELPER
function msg = make_failure_message(tflist, Tmin, Tmax)
   fails = nnz(~tflist);  % count the falses
   if fails == 1
      idx = find(fails, 1);
      msg = sprintf('Measurement #%d exeeds the specified temperature range [%g,%g]', idx, Tmin, Tmax);
   else
      msg = sprintf('%d measurements exceed the specified temperature range [%g,%g]', fails, Tmin, Tmax);
   end

end
