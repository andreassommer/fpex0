function fitsol = FPEX0_exampleFit()
% fitsol = FPEX0_exampleFit()
%
% Example for fitting with FPEX0 
%
% INPUT:  none
%
% OUTPUT: fitsol --> fitting solution object
%
%
% NOTES:
%   * make sure that the paths have been set by invoking FPEX0_initPaths()
%   * initialize a parallel computing pool, e.g. via parpool(8) if you have 8 CPU cores
%
%
% Andreas Sommer, Sep2022, Aug2026
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu
%

% load the example setup
FPEX0setup = FPEX0_exampleSetup();

% modify some configuration (as example)
%FPEX0setup.Integration.options.RelTol = 1.0d-8;
%FPEX0setup.Integration.options.AbsTol = 1.0d-14;
%FPEX0setup.Integration.VDEoptions.RelTol = 1.0d-8;
%FPEX0setup.Integration.VDEoptions.AbsTol = 1.0d-14;


% import the example data
FPEX0setup = FPEX0_importExampleMeasurements(FPEX0setup, 2); % 2 = gridskip


% solve the fitting problem
% lsqnonlin -> requires Matlab's Optimization Toolbox
% fminsearch -> poor derivative-free optimizer with manually forced non-negativity
% In FPEX0_fit.m, you can also add a substitute optimizer.
if which('lsqnonlin')
   optimizer = 'lsqnonlin';
else
   optimizer = 'fminsearch';
   warning("lsqnonlin() from Matlab's Optimization Toolbox not available. Fallback to %s", optimizer)
   input('Usually, this optimizer does not perform well for these type of problems. Press ENTER to try anyway.', 's');
end
tic
fitsol = FPEX0_fit(FPEX0setup, 'optimizer', optimizer);
toc

% finito
return

end % of function



