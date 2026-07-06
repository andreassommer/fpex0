function DSCsample = DSCtools_restrictToTrange(DSCsample, Tmin, Tmax)
% DSCsample = DSCtools_restrictToTrange(DSCsample, Tmin, Tmax)
%
% Restricts data in DSCsample to specified temperature range.
%
% INPUT:  DSCsample --> DSC data structure as returned by DSCtools_readFile
%              Tmin --> minimum temperature to restrict data to   [default: -inf]
%              Tmax --> maximum temperature to restrict data to   [default: +inf]
%
% OUTPUT: DSCsample --> DSC data structure with restricted temperature range
%
%
% Author: Andreas Sommer, Jun2026
% code@andreas-sommer.eu

% defaults
if (nargin < 2), Tmin = -inf(); end
if (nargin < 3), Tmax = +inf(); end

% recursive call if DSCsample is not scalar
if ~isscalar(DSCsample)
   DSCsample = arrayfun(@(x) DSCtools_restrictToTrange(x, Tmin, Tmax), DSCsample);
   return
end

% get the temperature indices and truncate the data
if ~all(isinf([Tmin Tmax]))
   T = DSCsample.data.T;
   % search the index from opposing sides, because we might have multiple times the same temperature in the beginning
   % due to residual from cooling rates or from "pumping"
   n = numel(T);
   if (T(1) == Tmin), idx_start = 1;  else,  idx_start = findFirstSmallerRev(T, Tmin, n, []); end
   if (T(n) == Tmax), idx_final = n;  else,  idx_final = findFirstGreater   (T, Tmax, 1, []); end
   % invalid temperature range specified?
   if isempty(idx_start) || isempty(idx_final)
      error('Specified temperature range [%g %g] exceeds data temperature range [%g %g]', Tmin, Tmax, min(T), max(T));
   end
   % truncate to specified range
   if (idx_start > 0) && (idx_final > idx_start) 
      DSCsample.data = truncate_arrays_in_struct(DSCsample.data, idx_start, idx_final);  % processes all subfields
   end
end

% finito
return

end % of function

