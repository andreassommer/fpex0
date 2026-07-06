function varargout = DSCtools_settings(preset)
% settings = DSCtools_settings(preset)
%
% Chooses a settings preset.
%
% INPUT:     preset --> named setting, e.g. 'default'. See code for details.
%
% OUTPUT:  settings --> settings structure
%
% GLOBAL SIDE EFFECT: Updates or replaces global DSCTOOLS_GLOBAL_SETTINGS
%
% Author:  Andreas Sommer, Aug2022, Jun2026
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu
%


% access to global settings structure
global DSCTOOLS_GLOBAL_SETTINGS         %#ok<GVMIS>


% choose setting
if (nargin == 0), error('No setting preset specified.'); end

% process command
switch upper(preset)
   case {'DEFAULT', 'DEFAULTS', 'ACTIVATE'}    % 'ACTIVATE' for backward compatibility
      settings = get_default_settings();

   case 'NETZSCH2026'
      settings = get_default_settings();
      settings.blds_preset = 'netzsch2026';
      % settings.addCP_scaleToRate = true;
      settings.addCP_scaleByMass = false;

   case 'RELOAD'
      preset = getSetting(DSCTOOLS_GLOBAL_SETTINGS, 'preset');
      DSCtools_settings(preset);
      return

   otherwise
      error('Unknown settings: %s', command)
end

% add name of preset
settings.preset = preset;

% set the global variable and return variable
DSCTOOLS_GLOBAL_SETTINGS = settings;
if (nargout > 0)
   varargout{1} = settings;
end

% finito
makeMessage('Activated global settings preset "%s"\n', preset);
return

end % of function

%% HELPERS








function defaults = get_default_settings()

% DEFAULT SETTINGS  ---  see the respective files for explanations
% ================

% addCP
defaults.addCP_debugMode    = true;
defaults.addCP_scaleByMass  = true;
defaults.addCP_scaleToRate  = true;
defaults.addCP_TrangeWarn   = false;
defaults.addCP_defaultSignalSource = 'uV';
defaults.addCP_defaultTmax  = 200;
defaults.addCP_defaultTmin  =  20;

% warnings if header columns (of data section) cannot be retrieved
defaults.dataMAT_warn_if_column_not_found  = true;
defaults.dataMAT_warn_if_column_not_unique = true;

% baseline settings 
defaults.getBaselinePrimitive_resolution = 1000; % resolution for baseline function (number of points)
defaults.getBaseline_haltOnWarning       = true; % halt if getBaseline has warnings

% baseline detection settings: blds -- choose preset
defaults.blds_preset = 'default';

% expected units
enc = DSCtools_getEncodingDefaults();  % get the encoding of the machine
defaults.expectedUnits = struct( ...
           'Time'        , 'min' , ...
           'Temperature' , sprintf('%cC', enc.deg), ...    % degree Celsius
           'MicroVolts'  , 'uV'  , ...
           'MilliWatts'  , 'mW'  , ...
           'Scaling'     , 'uV/mW' ...
           );

% for backup of settings
defaults.previousSettings = [];

end