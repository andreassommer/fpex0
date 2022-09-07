function FPEX0setup = FPEX0_exampleSetup()
% FPEX0setup = FPEX0_exampleSetup()
%
% Generates example setup for FPEX0.
%
% INPUT:  none
%
% OUTPUT: setup --> example setup (handle) object
%
% Author:  Andreas Sommer, 2017, 2018, Aug2022
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu
%


% NOTE:
% Diffusion and drift (their value) is connected to the grid step size!
% So simulations and fits should be done with identical gridT.
%
% Initial and final diffusion are described by parameters!  
% --> Linear diffusion plus condition of a non-negative diffusion at end time
%
% Formula: Diffusion = pD1 + t*(pD2-pD1)/betamax   
%
%   where: betamax --> maximum heat rate (maximum time in FP)
%          pD1     --> initial diffusion at heat rate (time) 0
%          pD2     --> final diffusion at heat rate (time) betamax


% NOTE:  Parameter bounds
% Parameter bounds should never be active in the solution. 
% Their only role is to limit search space and to ensure that no "invalid" valued are attained.


% Generate parameters object
p0_FPdrift       = [  0.1  ;  +0.1  ];    % linear drift
p_lb_FPdrift     = [ -0.1  ;  -0.1  ];    % lower bounds
p_ub_FPdrift     = [ 10.1  ;  10.1  ];    % upper bounds

p0_FPdiffusion   = [  0.2  ;   0.1  ];    % linear diffusion
p_lb_FPdiffusion = [  0.00 ;   0.00 ];    % lower bounds
p_ub_FPdiffusion = [ 10.00 ;  10.00 ];    % upper bounds

p0_iniDist       = [  2  ;  40  ; 135 ; 15.0 ; 0.1000 ];  % Fraser-Sizuki parameters
p_lb_iniDist     = [  1  ;  15  ; 110 ;  0.1 ; 1d-5   ];  % lower bounds
p_ub_iniDist     = [  3  ; 150  ; 150 ; 50.0 ; 1.00   ];  % upper bounds

parametersObj  = FPEX0_class_parameters(  p0_FPdrift,   p0_FPdiffusion,   p0_iniDist, ...
                                        p_lb_FPdrift, p_lb_FPdiffusion, p_lb_iniDist, ...
                                        p_ub_FPdrift, p_ub_FPdiffusion, p_ub_iniDist  );


% Generate grid 
N        = 501;  % space resolution
betamax  = 20;   % maximum heat rate
gridT    = linspace( 60,      160,  N );   % x-grid = temperatures
gridTdot = linspace(  0,  betamax, 20 );   % t-grid = heating rates
gridObj  = FPEX0_class_grid(gridT, gridTdot);


% set the FP functions and the initial distribution
FPdriftFcn     = @(t,p) FPEX0_defaultDriftFcn(t,p);
FPdiffusionFcn = @(t,p) FPEX0_defaultDiffusionFcn(t,p,betamax);
IniDistFcn     = @(x,p) frasersuzuki(x,p);


% Setup integrator (using defaults)
integrationObj = FPEX0_class_integration();


% generate the setup object
FPEX0setup = FPEX0_class_setup(gridObj, parametersObj, integrationObj, FPdriftFcn, FPdiffusionFcn, IniDistFcn);




%


end