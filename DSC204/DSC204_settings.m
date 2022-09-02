function settings = DSC204_settings(varargin)
% settings = DSC204_settings(varargin)
%
% Generates, updates, or replaces DSC204 settings.
%
%  INPUT:   key-value-pairs:
%         'replace' --> replaces global DSC204settings with new settings
%          'update' --> updates global DSC204settings with new settings
%
%  OUTPUT:  settings -> settings structure
%
% Special call:  DSC204_settings('ACTIVATE') loads and activates default settings
%
% Side effects: May update or replace global DSC204settings.
%
% Author:  Andreas Sommer, Aug2022
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu
%



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
% dataMAT layout
defaults.dataMAT_column_t  = 2;
defaults.dataMAT_column_T  = 1;
defaults.dataMAT_column_uV = 3;
defaults.dataMAT_column_sf = 4;

% =================
% Do not edit below



% access to global settings structure
global DSC204settings

% default values for input arguments
replaceCurrentSettings = false;
updateCurrentSettings = false;

% shortcut for "nice" call using single "ACTIVATE" keyword
if (nargin == 1) && strcmpi(varargin{1},'ACTIVATE')
   settings = DSC204_settings('replace',true);
   return
end

% process arguments
argList  = varargin;
argCount = length(argList);
if mod(argCount, 2) == 1; error('Invalid number of arguments'); end
% REPLACEMENT?
if hasOption(argList, 'replace')
   replaceCurrentSettings = getOption(argList, 'replace'); 
   argList  = removeOption(argList, 'replace');
   argCount = argCount - 2;
end
% UPDATES?
if hasOption(varargin, 'update')
   updateCurrentSettings = getOption(argList, 'update');
   argList  = removeOption(argList, 'update'); 
   argCount = argCount - 2;
end



% if settings are to be updated, use current settings instead of defaults
settings = defaults;
if updateCurrentSettings
   settings = DSC204settings;
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
   DSC204settings.previousSettings = [];        % remove old settings backup
   settings.previousSettings = DSC204settings;  % backup old settings
   DSC204settings = settings;                   % activate settings
end



