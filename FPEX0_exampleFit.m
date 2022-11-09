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
% Andreas Sommer, Sep2022
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
tic
fitsol = FPEX0_fit(FPEX0setup, 'optimizer', 'lsqnonlin');
toc

% finito
return

end



