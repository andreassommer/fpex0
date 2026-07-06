function DSCsample = DSCtools_compact(DSCsample, varargin)
% DSCsample = DSCtools_compact(DSCsample, key-value-pairs*)
%
% Removes unnecessary data from dsc data structure, e.g. raw data from reading input files.
% If specified, also a resampling to a unified temperature grid is done
%
% INPUT:  DSCsample --> DSC data structure as returned by DSCtools_readFile
%   key-value-pairs --> possible keys:
%                       Tmin: minimum temperature to restrict data to   [default: -inf]
%                       Tmax: maximum temperature to restrict data to   [default: +inf]
%                       dT:   resampling to temperature increment dT
%
% OUTPUT: DSCsample --> compacted DSC data structure
%
% NOTE:   If dT is empty or zero, no resampling will be done.
%
% Author: Andreas Sommer, Jun2026
% code@andreas-sommer.eu


% recursive call if DSCsample is not scalar
if ~isscalar(DSCsample)
   DSCsample = arrayfun(@(x) DSCtools_compact(x, varargin{:}), DSCsample);
   return
end

% check optional arguments
args = varargin;
[      Tmin, args] = olGetOption(args, 'Tmin'      , -inf() );
[      Tmax, args] = olGetOption(args, 'Tmax'      , +inf() );
[        dT, args] = olGetOption(args, 'dT'        , []     );
[cp_warning, args] = olGetOption(args, 'cp_warning', true   );
olWarnIfNotEmpty(args);

% delete not necessary fields
DSCsample.rawData.dataMat = [];    % raw data - we cannot call addQuickAccessors() anymore
DSCsample.rawData.fileEnc = [];    % file encoding - not needed

% restrict data to specified temperature
if ~all(isinf([Tmin Tmax]))
   T = DSCsample.data.T;
   idxT = isfinite(T) & (T >= Tmin) & (T <= Tmax);
   idx_start = find(idxT, 1, 'first');
   idx_final = find(idxT, 1, 'last');
   if (idx_start > 0) && (idx_final > idx_start) 
      DSCsample.data = truncate_arrays_in_struct(DSCsample.data, idx_start, idx_final);
   end
end

% resampling?
if ~isempty(dT)
   DSCsample = DSCtools_resample(DSCsample, dT, Tmin, Tmax, 'cp_warning', false);  % always suppress cp_warning
end

% if DSCsample contains cp field, issue a warning
if cp_warning
   if isfield(DSCsample, 'cp')
      makeMessage('cp field in DSCsample with ID "%s" is NOT UPDATED.\n')
   end
end

% finito
return
   
end % of function

