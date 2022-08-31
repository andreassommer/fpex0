function FPEX0_setup(varargin)
% FPEX0_setup(varargin)
%
% Creates FPEX0 configuration and stores it in global variable FPEX0
%
% INPUT:   key-value-pairs: none yet
%
% OUTPUT   config --> structure containing configuration
%
% Author:  Andreas Sommer, Aug2022
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu


% Remark: This would all fit a bit better inside an object oriented
%         framework using a class; however, that would limit compatibility.

% initialize setup variable
global FPEX0
FPEX0 = [];


% process arguments
args = varargin;
if hasOption(args, 'setup'), FPEX0 = getOption(args, 'setup'); end


% SET SOME DEFAULTS HERE
FPEX0.grid.N         = 1001;                                           % number of grid points
FPEX0.grid.gridT     = linspace(50, 160, FPEX0.grid.N);                % x-grid = temperatures
FPEX0.grid.gridTdot  = linspace( 0,  20, 20 );                         % t-grid = heating rates
FPEX0.grid.h         = diff(FPEX0.grid.gridT([end 1])) / FPEX0.grid.N; % grid step size
%setup.grid = grid;

% function handles
FPEX0.functions.drift      = @FPEX0_defaultDriftFcn;
FPEX0.functions.diffusion  = @FPEX0_defaultDiffusionFcn;
FPEX0.functions.make_rhs   = @make_FPODE_rhsfun;
FPEX0.functions.make_jac   = @make_FPODE_jacfun;
FPEX0.functions.jacPattern = @make_FPODE_jacpattern;
FPEX0.functions.iniDist    = @frasersuzuki;

% parameter layout
FPEX0.parameters.count              = 9;  % DEBUG  5IC + 2+2FP
FPEX0.parameters.values             = zeros(1, FPEX0.parameters.count);
FPEX0.parameters.idx.all            = 1:length(FPEX0.parameters.values);
FPEX0.parameters.idx.FPdrift        = 1:2;
FPEX0.parameters.idx.FPdiffusion    = 3:4;
FPEX0.parameters.idx.FPall          = [FPEX0.parameters.idx.FPdrift, FPEX0.parameters.idx.FPdiffusion]; % combined
FPEX0.parameters.idx.iniDist        = 5:9;

FPEX0.parameters.get.all            = @() getP_all();
FPEX0.parameters.get.FPall          = @() getP_FPall();
FPEX0.parameters.get.FPdrift        = @() getP_FPdrift();
FPEX0.parameters.get.FPdiffusion    = @() getP_FPdiffusion();
FPEX0.parameters.get.iniDist        = @() getP_iniDist();
FPEX0.parameters.set.all            = @setP_all;
FPEX0.parameters.set.FPall          = @setP_FPall;
FPEX0.parameters.set.FPdrift        = @setP_FPdrift;
FPEX0.parameters.set.FPdiffusion    = @setP_FPdiffusion;
FPEX0.parameters.set.iniDist        = @setP_iniDist;

FPEX0.parameters.extract.all         = @extractP_all;
FPEX0.parameters.extract.FPall       = @extractP_FPall;
FPEX0.parameters.extract.FPdrift     = @extractP_FPdrift;
FPEX0.parameters.extract.FPdiffusion = @extractP_FPdiffusion;
FPEX0.parameters.extract.iniDist     = @extractP_iniDist;


% This seemingly clumsy implementation ensures that we always access the current values,
% and not old ones over which has been closed
% Note: Think of this example:   FPEX0.parameters.get.all = @() FPEX0.parameters.values
%       Here, we close over FPEX0.parameters.values, i.e. these values are fixed in the instant
%       when the anonymous function @() is created! FPEX0 is never accessed afterwards!

% get internally stored parameters
function p = getP_all();         p = FPEX0.parameters.values;                                   end
function p = getP_FPall();       p = FPEX0.parameters.values(FPEX0.parameters.idx.FPall);       end
function p = getP_FPdrift();     p = FPEX0.parameters.values(FPEX0.parameters.idx.FPdrift);     end
function p = getP_FPdiffusion(); p = FPEX0.parameters.values(FPEX0.parameters.idx.FPdiffusion); end
function p = getP_iniDist();     p = FPEX0.parameters.values(FPEX0.parameters.idx.iniDist);     end

% set internally stored parameters
function setP_all(p);         FPEX0.parameters.values = p;                                   end
function setP_FPall(p);       FPEX0.parameters.values(FPEX0.parameters.idx.FPall) = p;       end
function setP_FPdrift(p);     FPEX0.parameters.values(FPEX0.parameters.idx.FPdrift) = p;     end
function setP_FPdiffusion(p); FPEX0.parameters.values(FPEX0.parameters.idx.FPdiffusion) = p; end
function setP_iniDist(p);     FPEX0.parameters.values(FPEX0.parameters.idx.iniDist) = p;     end

% extract parameters given as variable
function p = extractP_all(pp);         p = pp;                                   end
function p = extractP_FPall(pp);       p = pp(FPEX0.parameters.idx.FPall);       end
function p = extractP_FPdrift(pp);     p = pp(FPEX0.parameters.idx.FPdrift);     end
function p = extractP_FPdiffusion(pp); p = pp(FPEX0.parameters.idx.FPdiffusion); end
function p = extractP_iniDist(pp);     p = pp(FPEX0.parameters.idx.iniDist);     end


   
% simulation data
FPEX0.simulation.sol = NaN;  % will be filled by FPEX0_simulate


% integration setting
FPEX0.integration.reltol     = 1.0e-06;
FPEX0.integration.abstol     = 1.0e-10;
FPEX0.integration.integrator = @ode15s;
FPEX0.integration.updateOptionsJacobian = @update_FPODE_integratorOptions_Jacobian;
FPEX0.integration.options = make_FPODE_defaultIntegratorOptions(FPEX0);




end




%% HELPERS



% Jacobian dudf sparsity pattern
function jacpattern = make_FPODE_jacpattern(N)
   e = ones(N,1);
   jacpattern = spdiags([e e e], -1:1, N, N); 
end

function jacfun = make_FPODE_jacfun(h, driftFcn, driftParams, diffusionFcn, diffusionParams)
   jacfun = @(t,u) FokkerPlanckODE_dfdu(t, u, h, driftFcn, driftParams, diffusionFcn, diffusionParams);
end

function rhsfun = make_FPODE_rhsfun(h, driftFcn, driftParams, diffusionFcn, diffusionParams)
   rhsfun = @(t,u) FokkerPlanckODE(t, u, h, driftFcn, driftParams, diffusionFcn, diffusionParams);
end





% options and integrator selection % DEBUG / TODO:  must be adjusted when N changes
function opts = make_FPODE_defaultIntegratorOptions(FPEX0)
   opts = odeset(...
        'AbsTol'      , FPEX0.integration.abstol  ...
      , 'RelTol'      , FPEX0.integration.reltol  ...
      , 'BDF'         , 'off'    ...
      , 'vectorized'  , 'on'     ...
      , 'JPattern'    , make_FPODE_jacpattern(FPEX0.grid.N)  ...
      , 'NonNegative' , 1:FPEX0.grid.N            ...
      , 'stats'       , 'off'      );
end

% update options:  must be adapted when 
function opts = update_FPODE_integratorOptions_Jacobian(opts, jacobian)
   opts = odeset(opts, 'Jacobian', jacobian);
end

