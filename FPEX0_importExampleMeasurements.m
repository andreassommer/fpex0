function FPEX0setup = FPEX0_importExampleMeasurements(FPEX0setup, gridskip)
% FPEX0setup = test_FPEX0_importMeasurements()
% FPEX0setup = test_FPEX0_importMeasurements(FPEX0setup)
% FPEX0setup = test_FPEX0_importMeasurements(FPEX0setup, gridskip)
%
% Imports example data from file ExampleMeasurements.dat into an FPEX0setup object.
%
% INPUT:  FPEX0setup --> FPEX0 setup object to be used      [default: new from FPEX0_exampleSetup()]
%           gridskip --> use every n-th grid point for data [default: 1]
%
% OUTPUT: FPEX0setup --> example setup (handle) object with added example measurement data
%
% Andreas Sommer, Sep2022
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu
%

% use example setup if not specified
if (nargin < 1)
   FPEX0setup = FPEX0_exampleSetup();
end   

% set default gridskip if not specified
if (nargin < 2)
   gridskip = 1;  % take every n-th sample only
end


% load example data from file
fprintf('Loading example data . . . ');
exdatafile = 'ExampleMeasurements.mat';
exdata = load(exdatafile,'exdata'); % load data
exdata = exdata.exdata;             % unwrap data
fprintf('Loaded.\n');

% select only specific experiment by ID
targetID  = '16-407';
targetIdx = arrayfun(@(x) strcmp(x.ID,targetID), exdata);
exdata    = exdata(targetIdx);

% register data in FPEX0setup
fprintf('Importing example measurements to FPEX0setup . . .');
for k = 1:length(exdata)
   FPEX0setup.importMeasurements(exdata(k).ID, exdata(k).rate, exdata(k).cp.T, exdata(k).cp.latentdata, gridskip);
   % cp.latentdata contains the cp values without baseline
end
fprintf('Imported.\n');

