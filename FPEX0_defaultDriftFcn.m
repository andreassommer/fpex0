function fval = FPEX0_defaultDriftFcn(t,p)
% fval = FPEX0_defaultDriftFcn(t,p)
%
% Default drift function used in FPEX0
%
% INPUT:   t --> current time / heating rate
%          p --> drift parameter vector (drift-parameters only!)
% 
% OUTPUT   fval --> function value for drift
%
% Author:  Andreas Sommer, Aug2022
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu


% linear drift parametrization
fval = p(1) + p(2) * t;


end