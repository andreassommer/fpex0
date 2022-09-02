% Tests import of measurements

% uses example measuremet data


% access global data
global FPEX0

% load example data
fprintf('Loading example data . . . ');
exdatafile = fullfile(FPEX0.paths.base, 'ExampleMeasurements.mat');
exdata = load(exdatafile,'exdata'); % load data
exdata = exdata.exdata;             % unwrap data
fprintf('Loaded.\n');

% select only specific experiment by ID
targetID = '16-407';
targetIdx = arrayfun(@(x) strcmp(x.ID,'16-407'), exdata);
exdata = exdata(targetIdx);

% take every n-th sample only
skip = 1;

% Backup previous measurement data
if ~isempty(FPEX0.measurements)
   fprintf('Clearing existing measurement data . . .');
   measdatabackup = FPEX0.measurements;
   FPEX0.measurements = [];
   fprintf('Cleared.\n  -- Backup created in workspace variable "measdatabackup" \n');
end

% register data 
fprintf('Importing example measurements to FPEX0 . . .');
for k = 1:length(exdata)
   FPEX0_importMeasurements(exdata(k).ID, exdata(k).rate, exdata(k).cp.T, exdata(k).cp.latentdata, skip);
   % cp.latentdata contains the cp values without baseline
end
fprintf('Imported.\n');

