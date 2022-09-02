% Tests the FPEX0_calcresvec
%

% load example measurements
test_FPEX0_importMeasurements

% evaluate residual vector on nominal parameter values
fprintf('Evaluating residual vector . . . ');
params = FPEX0.parameters.get.all();
resvec = FPEX0_calcresvec(params);
fprintf('Evaluated.\n');
