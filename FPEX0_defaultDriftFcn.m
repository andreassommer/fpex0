function [fval, dfdp] = FPEX0_defaultDriftFcn(t,p)
% fval = FPEX0_defaultDriftFcn(t,p)
% [fval, dfdp] = FPEX0_defaultDriftFcn(t,p)
%
% Default drift function used in FPEX0
%
% INPUT:   t --> current time / heating rate
%          p --> drift parameter vector (drift-parameters only!)
% 
% OUTPUT   fval --> function value for drift
%          dfdp --> derivative w.r.t. p
%
% Author:  Andreas Sommer, Aug-Sep2022
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu


% linear drift parametrization
fval = p(1) + p(2) * t;

% derivative
if (nargout >= 2)
   t = reshape(t, 1, []);    % reshape such that dfdp(i,:)
   I = ones(size(t));        %
   dfdp = [ I ; t];          % contains derivatives for p_i
end


end