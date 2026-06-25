function settings = DSCtools_settings(varargin)
% settings = DSCtools_settings(varargin)
%
% Generates, updates, or replaces DSCtools settings.
%
%  INPUT:   key-value-pairs:
%         'replace' --> replaces global DSCTOOLS_GLOBAL_SETTINGS with new settings
%          'update' --> updates global DSCTOOLS_GLOBAL_SETTINGS with new settings
%
%  OUTPUT:  settings -> settings structure
%
% Special call:  DSCtools_settings('ACTIVATE') loads and activates default settings
%
% Side effects: May update or replace global DSCTOOLS_GLOBAL_SETTINGS.
%
% Author:  Andreas Sommer, Aug2022, Jun2026
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu
%


% get the encoding of the machine
enc = DSCtools_getEncodingDefaults();

% DEFAULT SETTINGS  ---  see the respective files for explanations
% ================

% addCP
defaults.addCP_debugMode   = true;
defaults.addCP_scaleByMass = true;
defaults.addCP_scaleToRate = true;
defaults.addCP_TrangeWarn  = false;
defaults.addCP_defaultSignalSource = 'uV';
defaults.addCP_defaultTmax = 160;
defaults.addCP_defaultTmin =  55;

% warnings if header columns (of data section) cannot be retrieved
defaults.dataMAT_warn_if_column_not_found  = true;
defaults.dataMAT_warn_if_column_not_unique = true;

% expected units
defaults.expectedUnits = struct( ...
           'Time'        , 'min' , ...
           'Temperature' , sprintf('%cC', enc.deg), ...    % degree Celsius
           'MicroVolts'  , 'uV'  , ...
           'MilliWatts'  , 'mW'  , ...
           'Scaling'     , 'uV/mW' ...
           );


% =================
% Do not edit below


% access to global settings structure
global DSCTOOLS_GLOBAL_SETTINGS         %#ok<GVMIS>

% default values for input arguments
replaceCurrentSettings = false;
updateCurrentSettings = false;

% shortcut for "nice" call using single "ACTIVATE" keyword
if (nargin == 1) && strcmpi(varargin{1},'ACTIVATE')
   settings = DSCtools_settings('replace',true);
   return
end

% process arguments
argList  = varargin;
argCount = length(argList);
if mod(argCount, 2) == 1; error('Invalid number of arguments'); end
% REPLACEMENT?
if olHasOption(argList, 'replace')
   replaceCurrentSettings = olGetOption(argList, 'replace'); 
   argList  = olRemoveOption(argList, 'replace');
   argCount = argCount - 2;
end
% UPDATES?
if olHasOption(varargin, 'update')
   updateCurrentSettings = olGetOption(argList, 'update');
   argList  = olRemoveOption(argList, 'update'); 
   argCount = argCount - 2;
end



% if settings are to be updated, use current settings instead of defaults
settings = defaults;
if updateCurrentSettings
   settings = DSCTOOLS_GLOBAL_SETTINGS;
end


% PROCESS ARGS
for k = 1:2:argCount
   argName = argList{k};    % get setting name
   argVal  = argList{k+1};  % always exists, as argCount is even (see check above)
   try % try to store the setting with specified name
      settings.(argName) = argVal;
   catch
      if     isnumeric(argName), error('Invalid setting at position #%g detected: "%g"', k, argName);
      elseif isstring(argName) , error('Invalid setting at position #%g detected: "%s"', k, argName);
      elseif ischar(argName)   , error('Invalid setting at position #%g detected: "%s"', k, argName);
      else                     , error('Invalid setting at position #%g detected: class %s', k, class(argName));
      end
   end
end



% check if settings should replaced or updated
if replaceCurrentSettings || updateCurrentSettings
   DSCTOOLS_GLOBAL_SETTINGS.previousSettings = [];        % remove old settings backup
   settings.previousSettings = DSCTOOLS_GLOBAL_SETTINGS;  % backup old settings
   DSCTOOLS_GLOBAL_SETTINGS = settings;                   % activate settings
end



