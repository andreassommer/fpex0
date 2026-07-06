function DSCsample = DSCtools_smoothdata(DSCsample, varargin)
% DSCsample = DSCtools_smoothdata(DSCsample, key-value-pairs)
%
% Applies a smoothing to the data in DSCsample.
% Temperature and time fields are not modified, and an equidistant stepping is assumed.
%
% INPUT:  DSCsample --> DSC data structure as returned by DSCtools_readFile
%   key-value-pairs --> arguments passed to Matlab's smoothdata() function
%
% OUTPUT: DSCsample --> DSC data structure with smoothed data
%
%
% Author: Andreas Sommer, Jun2026
% code@andreas-sommer.eu


% recursive call if DSCsample is not scalar
if ~isscalar(DSCsample)
   DSCsample = arrayfun(@(x) DSCtools_smoothdata(x, varargin{:}), DSCsample);
   return
end

% input arguments to be passed
args = varargin;

% quick accessor for data
data = DSCsample.data;

% list of to-be-processed subfields
fields = fieldnames(data);
fields(strcmp(fields,'T')) = [];  % remove the temperature field from the list
fields(strcmp(fields,'t')) = [];  % remove the        time field from the list

% smooth it
for i = 1:length(fields)          % process all fields (except T)
   f = fields{i};
   data.(f) = smoothdata(data.(f), args{:});  % apply smoothing
end
DSCsample.data = data;            % rewrite all to struct


% finito
return
   
end % of function
