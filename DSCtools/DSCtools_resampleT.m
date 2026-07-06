function DSCsample = DSCtools_resampleT(DSCsample, dT, Tmin, Tmax, varargin)
% DSCsample = DSCtools_resampleT(DSCsampleT, dT, Tmin, Tmax, key-value-pairs)
%
% Resamples recorded data to specified temperature increment and range.
%
% INPUT:  DSCsample --> DSC data structure as returned by DSCtools_readFile
%                dT --> temperature increment                     [default: 0.05]
%              Tmin --> minimum temperature to restrict data to   [default: -inf]
%              Tmax --> maximum temperature to restrict data to   [default: +inf]
%
% OUTPUT: DSCsample --> compacted DSC data structure
%
% NOTE:  The data is guaranteed to span the whole temperature range [Tmin, Tmax].
%        T(1) always equals Tmin.
%        T(end) may be larger than Tmax:  T(end) = min { T | T := Tmin + k*dt >= Tmax, k = 1,2,3 ... }
%
% Author: Andreas Sommer, Jun2026
% code@andreas-sommer.eu

% default args
if (nargin < 2),   dT = 0.05;   end
if (nargin < 3), Tmin = -inf(); end
if (nargin < 4), Tmax = +inf(); end

% check optional arguments
args = varargin;
[cp_warning, args] = olGetOption(args, 'cp_warning', true);
olWarnIfNotEmpty(args);

% signal error if temperature range exceeds available data
DSCtools_assertTrange(DSCsample, Tmin, Tmax, 'errorOnFail', true);

% recursive call if DSCsample is not scalar
if ~isscalar(DSCsample)
   DSCsample = arrayfun(@(x) DSCtools_resampleT(x, dT, Tmin, Tmax), DSCsample);
   return
end

% --- here DSCsample is a SINGLE (1-element array) struct

% restrict data to specified temperature
DSCsample = DSCtools_restrictToTrange(DSCsample, Tmin, Tmax);

% default bounds for interpolation, ensuring that [Tmin, Tmax] is covered
T_lo = Tmin;
T_hi = (Tmin + ceil((Tmax-Tmin)/dT)*dT);
if isinf(Tmin), Tmin = min(DSCsample.data.T); T_lo = Tmin; end  % on infinite T limits
if isinf(Tmax), Tmax = max(DSCsample.data.T); T_hi = Tmax; end  % we adjust the interpolation bounds
newT = T_lo:dT:T_hi;

% resample the fields
data = DSCsample.data;
fields = fieldnames(data);        % list of to-be-processed subfields
fields(strcmp(fields,'T')) = [];  % remove the temperature field from the list
oldT = data.T;                    % keep old temperatures for interpolation
data.T = newT;                    % store new temperature grid
for i = 1:length(fields)          % process all fields (except T)
   f = fields{i};
   data.(f) = interp1(oldT, data.(f), newT, 'linear');  % linear interpolation at new temperature grid
end
DSCsample.data = data;            % rewrite all to struct

% if DSCsample contains cp field, issue a warning
if cp_warning
   if isfield(DSCsample, 'cp')
      makeMessage('cp field in DSCsample with ID "%s" is NOT UPDATED.\n', DSCsample.ID);
   end
end

% finito
return
   
end % of function
