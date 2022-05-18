function [mfun, tvec] = DSC204test_makeMeasurementsFunctionFromMatFile(matfilename)
   % [mfun, tvec] = DSC204test_makeMeasurementsFunctionFromMatFile(matfilename)
   %
   % Creates measurement function from data in specified .mat file.
   %
   % INPUT:    matfilename --> file name of .mat file with measurements
   %                          (as created by DSC204test_makeMeasurementsMatFile.m)
   %
   % OUTPUT:          mfun --> function accepting time and spacial coordinates
   %                  tvec --> vector of available measurement times
   % 

   % load the mat file
   mdata = load(matfilename);
   mdata = mdata.mdata;  % resolve one level of nesting
   tvec  = mdata.tvec;   % measurement times
   
   T = mdata.T;               % T = temperature  (= spacial coordinate)
   lcpdata = mdata.latentcp;  % lcp = latent cp
   
   rates = mdata.tvec;   % time points = heating rates
   
   % nested function for measurement retrieval
   function eta = measf(t, x, epsilon)
      % t = rate, x = temperature
      if (nargin<3), epsilon = 0.01; end                      % default tolerance for finding time (rate)
      idx = find( (abs(rates - t) < epsilon), 1, 'first');    % index of data closest to specified time
      if isempty(idx)
         validtimes = sprintf(' %g ', rates);
         error('Invalid heating rate (FP time): %g.\nMust be one of [%s].\n', t, validtimes);
      end
      eta = interp1(T{idx}, lcpdata{idx}, x, 'linear', 0);
   end

   % output function
   mfun = @(t,x,epsilon) measf(t,x,epsilon);
   
end