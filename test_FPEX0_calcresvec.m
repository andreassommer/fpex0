function resvec = test_FPEX0_calcresvec()
% Tests the FPEX0_calcresvec

% load example measurements
gridskip = 1;
FPEX0setup = test_FPEX0_importMeasurements(gridskip);

% evaluate residual vector on nominal parameter values
fprintf('Evaluating residual vector . . . ');
p0     = FPEX0setup.Parameters.p0;
resvec = FPEX0_calcresvec(FPEX0setup, p0);
fprintf('Evaluated.\n');

% finito
return

end