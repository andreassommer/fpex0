function FPEX0setup = test_FPEX0_importMeasurements(gridskip)
% FPEX0setup = test_FPEX0_importMeasurements()
% FPEX0setup = test_FPEX0_importMeasurements(gridskip)
%
% Tests import of measurements using example data from file ExampleMeasurements.dat
%
% INPUT:  gridskip --> use every n-th grid point for data [default: 1]
%
% OUTPUT: FPEX0setup --> example setup (handle) object
%
% Author:  Andreas Sommer, 2017, 2018, Aug-Sep2022
% andreas.sommer@iwr.uni-heidelberg.de
% code@andreas-sommer.eu
%

% set default gridskip if not specified
if (nargin==0)
   gridskip = 1;  % take every n-th sample only
end

% generate example Setup
FPEX0setup = FPEX0_exampleSetup();

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

