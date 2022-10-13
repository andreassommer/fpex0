function [fval, dfdp] = FPEX0_defaultDiffusionFcn(t,p,betamax)
% fval = FPEX0_defaultDiffusionFcn(t,p,betamax)
% [fval, dfdp] = FPEX0_defaultDiffusionFcn(t,p,betamax)
%
% Default diffusion function used in FPEX0
%
% INPUT:      t --> current time / heating rate
%             p --> diffusion parameter vector (diffusion-parameters only!)
%       betamax --> maximum time / heat rate (used for ensuring non-negativity)
%                   (optional, default = 100);
% 
% OUTPUT   fval --> function value for diffusion
%          dfdp --> derivative w.r.t. p
%

% Author:  Andreas Sommer, Aug-Sep2022
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu


% linear parametrization ensuring non-negativity for non-negative p1 and p2
fval = p(1) + t * (p(2) - p(1)) / betamax;

% derivative
if (nargout >= 2)
   t = reshape(t, 1, []);                    % reshape such that dfdp(i,:) 
   dfdp = [ 1 - t/betamax ; t/betamax ];     % contains derivatives for p_i
end

end