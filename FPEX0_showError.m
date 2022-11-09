function FPEX0_showError(msg, err)
% FPEX0_showError(msg)
% FPEX0_showError(msg, err)
%
% Displays error message.
%
% INPUT:   msg --> message to be displayed
%          err --> matlab error structure to be decoded
%
% OUTPUT:  none
%
% Andreas Sommer, Nov2022
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu
%


fprintf('ERROR: %s\n', msg);
if (nargin >= 2)
   fprintf('Error id    : %s\n', err.identifier);
   fprintf('Error string: %s\n', err.message);
end


end