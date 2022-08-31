function fval = FPEX0_defaultDiffusionFcn(t,p,betamax)
% fval = FPEX0_defaultDiffusionFcn(t,p,betamax)
%
% Default diffusion function used in FPEX0
%
% INPUT:      t --> current time / heating rate
%             p --> diffusion parameter vector (diffusion-parameters only!)
%       betamax --> maximum time / heat rate (used for ensuring non-negativity)
%                   (optional, default = 100);
% 
% OUTPUT   fval --> function value for diffusion
%
% Author:  Andreas Sommer, Aug2022
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu

% use default maxheatrate if not specified
if (nargin<3)
   betamax = 100;
end

% linear parametrization ensuring non-negativity for non-negative p1 and p2
fval = p(1) + t * (p(2) - p(1)) / betamax;



end