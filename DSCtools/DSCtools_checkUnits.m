function allOK = DSCtools_checkUnits(varargin)
% allOK = DSCtools_checkUnits(key-value-pairs*)
%
% Verifies that the units are as expected.
%
% INPUT:  key-value-pairs with possible keys:
%         'Time'        --> variable carrying the unit of time
%         'Temperature' --> variable carrying the unit of temperature
%         'MicroVolts'  --> variable carrying the unit of micro Volts
%         'MilliWatts'  --> variable carrying the unit of milli Watts
%         'Scaling'     --> variable carrying the unit of scaling factor, i.e. mW/uV
% 
% OUTPUT: allOK --> true if all units are as expected, false otherwise
%
% NOTE:  If an input is empty [], it is not checked.
%        If an input is empty '', is will be considered as failed.
%
% Author:   Jun2026
% code@andreas-sommer.eu

% retrieve expected units from settings
global DSCTOOLS_GLOBAL_SETTINGS  %#ok<GVMIS>
expectedUnits = getSetting(DSCTOOLS_GLOBAL_SETTINGS, 'expectedUnits');


% we return as soon as a test failed
allOK = false;

% get args
args = varargin;
[ t, args] = olGetOption(args, 'time'       , []);
[ T, args] = olGetOption(args, 'Temperature', []);
[uV, args] = olGetOption(args, 'MicroVolts' , []);
[mW, args] = olGetOption(args, 'MilliWatts' , []);
[sf, args] = olGetOption(args, 'Scaling'    , []);
olWarnIfNotEmpty(args, true); % something left?

if ~sameUnits(  t, expectedUnits.Time        ), return, end
if ~sameUnits(  T, expectedUnits.Temperature ), return, end
if ~sameUnits( uV, expectedUnits.MicroVolts  ), return, end
if ~sameUnits( mW, expectedUnits.MilliWatts  ), return, end
if ~sameUnits( sf, expectedUnits.Scaling     ), return, end


% all tests passed, so okay
allOK = true;

% finito
return
   
end % of function



%% HELPERS
function tf = isSpecified(x)
   tf = ~ ( isnumeric(x) && isempty(x) );
end

function tf = sameUnits(x, expectedUnit)
   tf = true;
   if isSpecified(x)
      tf = strcmp(x, expectedUnit);
   end
   if ~tf
      makeMessage('Unit is not as expected: "%s" is not "%s\n', x, expectedUnit)
      keyboard
   end
end