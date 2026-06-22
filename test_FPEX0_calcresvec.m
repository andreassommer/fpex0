%function resvec = test_FPEX0_calcresvec()
% Tests the FPEX0_calcresvec

% load example measurements
FPEX0setup = FPEX0_exampleSetup();
FPEX0setup = FPEX0_importExampleMeasurements(FPEX0setup, 2); % 2 = gridskip


% evaluate residual vector on nominal parameter values
makeMessage('Evaluating residual vector . . . ');
p0     = FPEX0setup.Parameters.p0;
resvec = FPEX0_calcresvec(FPEX0setup, p0);
makeMessage('Evaluated.\n');
makeMessage('Norm of residual vector: %g\n', norm(resvec));
makeMessage('\n')

% evaluate with jacobian
makeMessage('Evaluating residual vector with jacobian. . . ');
p0     = FPEX0setup.Parameters.p0;
[resvec, jacobian] = FPEX0_calcresvec(FPEX0setup, p0);
makeMessage('Evaluated.\n');
makeMessage('Norm of residual vector: %g\n', norm(resvec));
makeMessage('\n')


% finito
return

%end